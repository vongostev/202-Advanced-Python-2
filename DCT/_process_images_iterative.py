# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 22:04:56 2021

@author: vonGostev
"""

import sys
import numpy as np
from scipy.linalg import lstsq
from loguru import logger as log
from time import perf_counter

__all__ = ('deconstruct_image', 'construct_image',
           'dct_matrix', 'filter_img_dct')


def deconstruct_image(
        img: np.ndarray,
        base_vecs: object) -> np.ndarray:
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
    img : np.ndarray
        Массив изображения.
    base_vecs : generator
        Набор базисных векторов (изображений).

    Returns
    -------
    base_coeffs : np.ndarray
        Набор коэффициентов разложения
        `img` по `base_vecs`.

    """
    img_flatten = img.flatten()
    t = perf_counter()
    base_coeffs = []
    for vec in base_vecs:
        base_coeffs.append(vec @ img_flatten)
    log.info(
        f'Image with shape {img.shape} deconstructed with matrix multiplication')
    log.info(f'Elapsed time {perf_counter() - t: .2f} s')
    return np.array(base_coeffs)


def construct_image(
        coeffs: np.ndarray,
        base_vecs: object,
        img_shape: tuple) -> np.ndarray:
    """
    Генерация изображения из набора коэффициентов разложения
    по набору базисных векторов.

    Parameters
    ----------
    coeffs : np.ndarray
        Набор коэффициентов разложения.
    base_vecs : object, generator
        Генератор матрицы базисных векторов преобразования.
    img_shape : tuple
        Размерность изображения.

    Returns
    -------
    img : np.ndarray
        Сгенерированное изображение.

    """
    img = sum(vec * coeff for vec, coeff in zip(base_vecs, coeffs))
    return np.array(img).reshape(img_shape)


def dct_matrix(dims: tuple) -> np.ndarray:
    """
    Генератор матрицы базисных векторов дискретного косинусного преобразования.

    Parameters
    ----------
    dims : tuple
        Размерность матрицы базисных векторов
        дискретного косинусного преобразования.

    Returns
    -------
    yield generator
        Генератор матрицы базисных векторов дискретного косинусного преобразования.

    """
    K, N = dims
    if K == -1:
        K = N
    n = np.arange(N)
    for k in np.arange(K).reshape((-1, 1)):
        if k > 0:
            yield np.sqrt(2 / N) * np.cos((np.pi * k * (1/2 + n)) / N)
        else:
            yield np.sqrt(1/N) * np.ones(N)


def filter_img_dct(
        img: np.ndarray, dct_dim: int = -1,
        filter_ampl: float = 10.):
    """
    Итеративная фильтрация пространственных частот изображения в базисе
    дискретного косинусного преобразования для больших изображений.
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

    Returns
    -------
    img_filtered : np.ndarray
        Изображение с фильтрованными спектральными компонентами.

    """
    img_flatten_shape = np.prod(img.shape)
    base_vecs = dct_matrix((dct_dim, img_flatten_shape))
    spectrum = deconstruct_image(img, base_vecs)
    spectrum[np.abs(spectrum) < filter_ampl] = 0
    log.info(f'DCT spectrum filtered by amplitudes less than {filter_ampl}')
    base_vecs = dct_matrix((dct_dim, img_flatten_shape))
    img_filtered = construct_image(spectrum, base_vecs, img.shape)
    log.info(
        f"Filtering made STD error {np.linalg.norm(img - img_filtered):.5g}")
    return img_filtered
