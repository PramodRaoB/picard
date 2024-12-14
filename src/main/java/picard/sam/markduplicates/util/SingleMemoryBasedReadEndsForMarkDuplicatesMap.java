/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package picard.sam.markduplicates.util;

import java.util.*;

/**
 * Map from String to ReadEnds object.  Memory-based implementation.  Used for MarkDuplicates.
 *
 * @author alecw@broadinstitute.org
 */
public class SingleMemoryBasedReadEndsForMarkDuplicatesMap implements ReadEndsForMarkDuplicatesMap {

    /**
     * Stores key, value pairs for a given sequence
     */
    private final Map<String, ReadEndsForMarkDuplicates> sequenceMap = new HashMap<>();

    public ReadEndsForMarkDuplicates remove(int mateSequenceIndex, String key) {
        return sequenceMap.remove(key);
    }

    public void put(int mateSequenceIndex, String key, ReadEndsForMarkDuplicates readEnds) {
        sequenceMap.put(key, readEnds);
    }

    public int size() {
        return sequenceMap.size();
    }

    /**
     * @return number of elements stored in RAM.  Always <= size()
     */
    public int sizeInRam() {
        return size();
    }

    public Set<Map.Entry<String, ReadEndsForMarkDuplicates>> entrySet() {
        return sequenceMap.entrySet();
    }
}
