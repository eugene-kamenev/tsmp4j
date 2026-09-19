package com.github.eugene.kamenev.tsmp4j.algo.extras.windowfinder;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class WindowFinderTest extends BaseTest {

    @Test
    void testWindowFinder() {
        var ts = data.stream()
            .mapToDouble(t -> t.x())
            .toArray();
        var result = MWF.mwf(ts, 10, 500);

        assertEquals(45, result);
    }
}
