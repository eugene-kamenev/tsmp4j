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

package com.github.eugene.kamenev.tsmp4j.algo.mp;

public class BaseMatrixProfile implements MatrixProfile {

    private final double[] profile;
    private final int[] indexes;

    private final double[] rightProfile;

    private final double[] leftProfile;

    private final int[] rightIndexes;

    private final int[] leftIndexes;

    public BaseMatrixProfile(double[] profile, int[] indexes) {
        this(profile, indexes, null, null, null, null);
    }

    public BaseMatrixProfile(double[] profile, int[] indexes, double[] rightProfile,
        double[] leftProfile, int[] rightIndexes, int[] leftIndexes) {
        this.profile = profile;
        this.indexes = indexes;
        this.rightProfile = rightProfile;
        this.leftProfile = leftProfile;
        this.rightIndexes = rightIndexes;
        this.leftIndexes = leftIndexes;
    }

    @Override
    public double[] profile() {
        return this.profile;
    }

    @Override
    public int[] indexes() {
        return this.indexes;
    }

    @Override
    public double[] rightProfile() {
        return this.rightProfile;
    }

    @Override
    public double[] leftProfile() {
        return this.leftProfile;
    }

    @Override
    public int[] rightIndexes() {
        return this.rightIndexes;
    }

    @Override
    public int[] leftIndexes() {
        return this.leftIndexes;
    }
}
