# -*- coding: utf-8 -*-
"""
Created on Sat Nov 27 21:08:46 2021

@author: vonGostev
"""
import numpy as np
import unittest
from cv2 import imread
from scipy.fft import dct
import _process_images as test
import _process_images_iterative as test_iter
import _process_images as good
import _process_images_iterative as good_iter


class DCTImageFilteringTest(unittest.TestCase):

    def setUp(self):
        self.images = [
            'img/logo64.png',
            'img/teapot128.png',
            'img/pcf_microphoto.jpg'
        ]
        self.test = test
        self.good = good

    def __test_spectrum_crop(self, path, dim, ampl=0, prec=1e-8):
        img = imread(path, 0)
        img_filtered = self.test.filter_img_dct(img, dim, filter_ampl=ampl)
        img_filtered_good = self.good.filter_img_dct(
            img, dim, filter_ampl=ampl)
        self.assertTrue(
            np.allclose(img_filtered, img_filtered_good, rtol=prec))

    def __test_deconstruct(self, path, dim, prec=1e-8):
        img = imread(path, 0)
        img_flatten_shape = np.prod(img.shape)
        bv = self.test.dct_matrix((dim, img_flatten_shape))
        coeffs = self.test.deconstruct_image(img, bv)
        bv = self.good.dct_matrix((dim, img_flatten_shape))
        coeffs_good = self.good.deconstruct_image(img, bv)
        self.assertTrue(
            np.allclose(coeffs, coeffs_good, rtol=prec))

    def test_dct_matrix(self):
        dims = np.arange(8, 1025, 8)
        for dim in dims:
            with self.subTest(i=dim):
                np.allclose(dct(np.eye(dim)),
                            np.array(list(self.test.dct_matrix((-1, dim)))).T)

    def test_deconstruct(self):
        modsnum = [1, 1000, 2000, 3000, 4000]
        for modnum in modsnum:
            for path in self.images:
                with self.subTest(i=modnum):
                    self.__test_deconstruct(path, modnum)

    def test_full_spectrum(self):
        for path in self.images:
            with self.subTest(i=path):
                img = imread(path, 0)
                img_filtered = self.test.filter_img_dct(
                    img, dct_dim=-1, filter_ampl=0)
                self.assertTrue(np.allclose(img_filtered, img))

    def test_filtering(self):
        modsnum = [1, 1000, 2000, 3000, 4000]
        ampls = [0, 1, 10, 100]
        for modnum in modsnum:
            for ampl in ampls:
                for path in self.images:
                    with self.subTest(i=f'{path}: Mods={modnum} AF={ampl}'):
                        self.__test_spectrum_crop(path, modnum, ampl)


class DCTIterativeImageFilteringTest(DCTImageFilteringTest):

    def setUp(self):
        super().setUp()
        self.test = test_iter
        self.good = good_iter


if __name__ == "__main__":
    unittest.main()
