/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.github.eugene.kamenev.tsmp4j;

import java.util.Objects;

/**
 * Mutable complex number class
 */
public class Complex {

    private double re;
    private double im;

    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }

    public double getReal() {
        return re;
    }
    public double getImaginary() {
        return im;
    }

    /**
     * @param c input
     * @return sum of this and c
     */
    public Complex add(Complex c) {
        re = re + c.re;
        im = im + c.im;
        return this;
    }

    /**
     * @param c input
     * @return multiplication of this and c
     */
    public Complex mult(Complex c) {
        Complex a = this;
        double real = a.re * c.re - a.im * c.im;
        double imag = a.re * c.im + a.im * c.re;
        return new Complex(real, imag);
    }

    @Override
    public String toString() {
        return re + ((Math.signum(im) >= 0) ? " + " : " - ") + Math.abs(im) + "i";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Complex complex = (Complex) o;
        return Double.compare(complex.re, re) == 0 && Double.compare(complex.im, im) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(re, im);
    }
}