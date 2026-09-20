package com.github.eugene.kamenev.tsmp4j.algo.pmp;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.utils.Util;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SKIMPTest extends BaseTest {

    @Test
    void testSkimpNoCrossCorrelation() {
        var limit = 200;
        var windows = Util.createRange(4, 6, 1);
        var skimp = new SKIMP(limit, false, windows);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(skimp::update);
        var result = skimp.get();
        var dst = result.profile();
        var ids = result.indexes();

        assertEquals(3, dst.length);
        assertEquals(3, ids.length);
        assertEquals(limit - 3, dst[0].length);
        assertEquals(limit - 4, dst[1].length);
        assertEquals(limit - 5, dst[2].length);
        assertEquals(limit - 3, ids[0].length);
        assertEquals(limit - 4, ids[1].length);
        assertEquals(limit - 5, ids[2].length);
    }

    @Test
    void testSkimpCrossCorrelation() {
        var limit = 200;
        var windows = Util.createRange(4, 6, 1);
        var skimp = new SKIMP(limit, true, windows);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(skimp::update);
        var result = skimp.get();
        var dst = result.profile();
        var ids = result.indexes();

        assertEquals(3, dst.length);
        assertEquals(3, ids.length);
        assertEquals(limit - 3, dst[0].length);
        assertEquals(limit - 4, dst[1].length);
        assertEquals(limit - 5, dst[2].length);
        assertEquals(limit - 3, ids[0].length);
        assertEquals(limit - 4, ids[1].length);
        assertEquals(limit - 5, ids[2].length);
    }
}
