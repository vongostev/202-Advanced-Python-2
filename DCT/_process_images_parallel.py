# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 22:04:56 2021

@author: vonGostev
"""

import sys
import numpy as np
from loguru import logger as log
from time import perf_counter
from joblib import Parallel, delayed
from cv2 import imread
from numba import njit, prange


@njit('float64[:, :](UniTuple(int64, 2), int64)',
      cache=True, nogil=True, fastmath=True, parallel=False)
def dct_matrix(dims: tuple, k_offset: int) -> np.ndarray:
    """
    Генерация матрицы базисных векторов дискретного косинусного преобразования.

    Parameters
    ----------
    dims : tuple
        Размерность матрицы базисных векторов
        дискретного косинусного преобразования.
    k_offset: int
        Номер начального базисного вектора.
        По умолчанию равен 0, то есть отчет векторов начинается с 0.
    Returns
    -------
    _dct_matrix : np.ndarray
        Матрица базисных векторов дискретного косинусного преобразования.

    """
    K, N = dims
    n = np.arange(N).reshape((1, N))
    matrix = np.zeros((N, K))
    for k in np.arange(K):
        if k + k_offset > 0:
            matrix[:, k] = np.sqrt(2 / N) * \
                np.cos((np.pi * (k + k_offset) * (1/2 + n)) / N)
        else:
            matrix[:, k] = np.sqrt(1 / N)
    return matrix


# @njit('float64[:](float64[:, :], float64[:])', fastmath=True, cache=True)
def _direct(base_vecs: np.ndarray, img_flatten: np.ndarray) -> np.ndarray:
    """
    Генерация коэффициентов разложения изображения по базисным векторам.

    Parameters
    ----------
    base_vecs : np.ndarray
        Матрица базисных векторов преобразования.
    img_flatten : np.ndarray
        Изображение в форме 1D массива.

    Returns
    -------
    C : np.ndarray
        Коэффициенты разложения.

    """
    return base_vecs.T.dot(img_flatten)


# @njit('float64[:](float64[:, :], float64[:])', cache=True, fastmath=True)
def _inverse(base_vecs: np.ndarray, spectrum: np.ndarray) -> np.ndarray:
    """
    Восстановление изображения по коэффициентам разложения

    Parameters
    ----------
    base_vecs : np.ndarray
        Матрица базисных векторов преобразования.
    spectrum : np.ndarray
        Коэффициенты разложения.

    Returns
    -------
    img : np.ndarray
        Изображение в форме 1D массива.

    """
    return base_vecs.dot(spectrum)


class DCTFilter:

    def __init__(self, path: str,
                 dct_dim: int = -1,
                 chunk_size: int = -1,
                 n_jobs: int = 1):
        """
        Конструктор класса фильтрации изображений с помощью
        дискретного косинусного преобразования (ДКП)

        Parameters
        ----------
        path : str
            Путь к изображению.
        dct_dim : int, optional
            Число используемых базисных векторов ДКП.
            По умолчанию -1, то есть используются все.
        chunk_size : int, optional
            Размер чанка в параллельных вычислениях.
            По умолчанию -1, то есть вычисляется сразу все.
        n_jobs : int, optional
            Число процессов параллельных вычислений. По умолчанию 1.

        Raises
        ------
        ValueError
            Путь к изображению некорректен.

        """

        self.path = path
        self.n_jobs = n_jobs

        self.img = imread(self.path, 0)
        if self.img is None:
            raise ValueError(f'Path `{path}` does not contain correct image')
        self.img_flatten = self.img.flatten().astype(np.float64)
        self.img_flatten_shape = len(self.img_flatten)

        if dct_dim == -1:
            dct_dim = self.img_flatten_shape
        if chunk_size == -1:
            offsets = [0]
            chunk_size = dct_dim
        else:
            offsets = np.arange(0, dct_dim, chunk_size)
        self.chunk_size = chunk_size
        self.dct_dim = dct_dim
        self.offsets = offsets

    def deconstruct(
            self, base_vecs: np.ndarray) -> np.ndarray:
        """
        Генерация коэффициентов разложения изображения
        в градациях серого по базисным векторам
        некоторого преобразования.

        Для минимизации отклонения изображения
        в базисе некоторого преобразования относительно
        исходного используется матричный метод
        для поиска коэффициентов разложения.

        Parameters
        ----------
        base_vecs : np.ndarray
            Набор базисных векторов (изображений).

        Returns
        -------
        base_coeffs : np.ndarray
            Набор коэффициентов разложения
            `img` по `base_vecs`.

        """
        t = perf_counter()
        base_coeffs = _direct(base_vecs, self.img_flatten)
        log.info(
            f'Image with shape {self.img.shape} deconstructed with matrix multiplication')
        log.info(f'Elapsed time {perf_counter() - t: .2f} s')
        return base_coeffs

    def construct(
            self,
            base_vecs: np.ndarray, k_offset: int = 0) -> np.ndarray:
        """
        Генерация изображения из набора коэффициентов разложения
        по набору базисных векторов.

        Parameters
        ----------
        base_vecs : np.ndarray
            Матрица базисных векторов преобразования.
        k_offset: int
            Номер начального базисного вектора.
            По умолчанию равен 0, то есть отчет векторов начинается с 0.

        Returns
        -------
        img : np.ndarray
            Сгенерированное изображение.

        """
        chunk_size = base_vecs.shape[1]
        img = _inverse(
            base_vecs, self.spectrum[k_offset:k_offset + chunk_size])
        return img

    def filter(self, filter_ampl: float = 0):
        """
        Итеративная фильтрация пространственных частот изображения в базисе
        дискретного косинусного преобразования для больших изображений.

        Parameters
        ----------
        filter_ampl : float, optional
            Амплитуда фильтрации. The default is 0.

        """

        def __d(k_offset):
            chunk_size = min(self.chunk_size, self.dct_dim - k_offset)
            base_vecs = dct_matrix(
                (chunk_size, self.img_flatten_shape), k_offset)
            return self.deconstruct(base_vecs)

        t = perf_counter()
        self.spectrum = np.hstack(
            Parallel(n_jobs=self.n_jobs,
                     prefer='processes')(
                delayed(__d)(x) for x in self.offsets))
        log.info(f'Deconstruction. Elapsed time {perf_counter() - t} s')
        self.spectrum[np.abs(self.spectrum) < filter_ampl] = 0
        log.info(
            f'DCT spectrum filtered by amplitudes less than {filter_ampl}')

        def __c(k_offset):
            chunk_size = min(self.chunk_size, self.dct_dim - k_offset)
            base_vecs = dct_matrix(
                (chunk_size, self.img_flatten_shape), k_offset)
            return self.construct(base_vecs, k_offset)

        t = perf_counter()
        self.__imgf = np.sum(
            Parallel(n_jobs=self.n_jobs,
                     prefer='processes')(
                         delayed(__c)(x) for x in self.offsets), axis=0)
        log.info(f'Construction. Elapsed time {perf_counter() - t} s')

    @property
    def filtered(self):
        """
        Возврат сгенерированного фильтрованного изображения

        Returns
        -------
        __imgf: np.ndarray
            Сгенерированное изображение.

        """
        return self.__imgf.reshape(self.img.shape)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    dctf = DCTFilter('img/pcf_microphoto.jpg', chunk_size=1024, n_jobs=4)
    dctf.filter()

    fig, axes = plt.subplots(1, 2)
    axes[0].imshow(dctf.img)
    axes[1].imshow(dctf.filtered)
    plt.show()
