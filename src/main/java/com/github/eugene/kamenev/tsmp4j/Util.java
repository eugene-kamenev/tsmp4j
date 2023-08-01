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

public class Util {
    public static double[] reverse(double[] arr) {
        double[] inv = new double[arr.length];
        for (int i=0; i<inv.length; i++) {
            inv[i] = arr[arr.length-1-i];
        }
        return inv;
    }

    public static int[] createRange(int start, int end, int step) {
        // Create the range of integers in the specified interval and step
        int length = (end - start) / step + 1;
        int[] range = new int[length];
        for (int i = 0; i < length; i++) {
            range[i] = start + i * step;
        }
        return range;
    }
}
