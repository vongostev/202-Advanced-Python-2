# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:44:31 2021

@author: von.gostev
"""
import sys
import numpy as np
from scipy.linalg import lstsq
from logging import Logger, StreamHandler, Formatter
from time import perf_counter

__all__ = ('deconstruct_image', 'construct_image',
           'dct_matrix', 'filter_img_dct')

log = Logger('dct.imgs')

# handler = StreamHandler(sys.stderr)
# handler.setLevel(30)
# formatter = Formatter(
#     '%(asctime)s - %(name)-10.10s [%(levelname)-7.7s]  %(message)s')
# handler.setFormatter(formatter)
# log.addHandler(handler)


def deconstruct_image(
        img: np.ndarray,
        base_vecs: np.ndarray,
        method: str = 'mul') -> np.ndarray:
    """
    Генерация коэффициентов разложения изображения
    в градациях серого по базисным векторам
    некоторого преобразования.

    Для минимизации отклонения изображения
    в базисе некоторого преобразования относительно
    исходного используется матричный метод или метод наименьших квадратов
    для поиска коэффициентов разложения.

    Parameters
    ----------
    img : np.ndarray
        Массив изображения.
    base_vecs : np.ndarray
        Набор базисных векторов (изображений).

    Returns
    -------
    base_coeffs : np.ndarray
        Набор коэффициентов разложения
        `img` по `base_vecs`.

    """
    assert method in ['mul', 'lstsq']

    img_flatten = img.flatten()
    t = perf_counter()
    if method == 'mul':
        base_coeffs = base_vecs.T @ img_flatten
        log.info(
            f'Image with shape {img.shape} deconstructed with matrix multiplication')
    elif method == 'lstsq':
        # Возвращаем только коэффициенты разложения
        base_coeffs = lstsq(base_vecs, img_flatten)[0]
        log.info('Image deconstructed with LSTSQ')
    log.info(f'Elapsed time {perf_counter() - t: .2f} s')
    return base_coeffs


def construct_image(
        coeffs: np.ndarray,
        base_vecs: np.ndarray,
        img_shape: tuple) -> np.ndarray:
    """
    Генерация изображения из набора коэффициентов разложения
    по набору базисных векторов.

    Parameters
    ----------
    coeffs : np.ndarray
        Набор коэффициентов разложения.
    base_vecs : np.ndarray
        Матрица базисных векторов преобразования.
    img_shape : tuple
        Размерность изображения.

    Returns
    -------
    img : np.ndarray
        Сгенерированное изображение.

    """
    img = base_vecs.dot(coeffs).reshape(img_shape)
    return img


def dct_matrix(dims: tuple) -> np.ndarray:
    """
    Генерация матрицы базисных векторов дискретного косинусного преобразования.

    Parameters
    ----------
    dims : tuple
        Размерность матрицы базисных векторов
        дискретного косинусного преобразования.

    Returns
    -------
    _dct_matrix : np.ndarray
        Матрица базисных векторов дискретного косинусного преобразования.

    """
    K, N = dims
    if K == -1:
        K = N
    k = np.arange(K).reshape((-1, 1))
    n = np.arange(N)
    _dct_matrix = np.sqrt(2 / N) * np.cos((np.pi * k * (1/2 + n)) / N)
    _dct_matrix[0, :] = np.sqrt(1/N)
    return _dct_matrix.T


def filter_img_dct(
        img: np.ndarray, dct_dim: int = -1,
        filter_ampl: float = 10., method: str = 'mul'):
    """
    Фильтрация пространственных частот изображения в базисе
    дискретного косинусного преобразования.
    Имитация JPEG сжатия

    Parameters
    ----------
    img : np.ndarray
        Массив изображения.
    dct_dim : int, optional
        Число векторов для разложения. The default is -1.
    filter_ampl : float, optional
        Амплитуда фильтрации.
        Если амплитуда компоненты меньше нее, то она зануляется.
        The default is 10.
    method : str, optional = ['mul', 'lstsq']
        Метод разложения изображения. The default is 'mul'.

    Returns
    -------
    img_filtered : np.ndarray
        Изображение с фильтрованными спектральными компонентами.

    """
    img_flatten_shape = np.prod(img.shape)
    t = perf_counter()
    base_vecs = dct_matrix((dct_dim, img_flatten_shape))
    log.info(
        f'{dct_dim} DCT vectors generated. Elapsed time {perf_counter() - t:.2f} s')
    spectrum = deconstruct_image(img, base_vecs, method)
    return spectrum
    spectrum[np.abs(spectrum) < filter_ampl] = 0
    log.info(f'DCT spectrum filtered by amplitudes less than {filter_ampl}')
    img_filtered = construct_image(spectrum, base_vecs, img.shape)
    log.info(
        f"Filtering made STD error {np.linalg.norm(img - img_filtered):.5g}")
    return img_filtered
