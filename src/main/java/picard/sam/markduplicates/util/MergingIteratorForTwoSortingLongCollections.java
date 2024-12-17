package picard.sam.markduplicates.util;

import htsjdk.samtools.util.*;

import java.util.*;

/**
 * Very specific utility class that can merge at most two sorting long collections
 */
public class MergingIteratorForTwoSortingLongCollections {
    SortingLongCollectionForMarkDuplicates first, second;
    boolean firstDone, secondDone;
    private static final Log log = Log.getInstance(MergingIteratorForTwoSortingLongCollections.class);

    public MergingIteratorForTwoSortingLongCollections(SortingLongCollectionForMarkDuplicates first, SortingLongCollectionForMarkDuplicates second) {
        this.first = first;
        this.second = second;
        if (this.first.hasNext()) firstDone = false;
        else {
            firstDone = true;
            first.cleanup();
        }
        if (this.second.hasNext()) secondDone = false;
        else {
            secondDone = true;
            second.cleanup();
        }
    }

    public MergingIteratorForTwoSortingLongCollections(SortingLongCollectionForMarkDuplicates first) {
        this.first = first;
        this.second = null;
        if (this.first.hasNext()) firstDone = false;
        else {
            firstDone = true;
            first.cleanup();
        }
        secondDone = true;
    }

    public boolean hasNext() {return !firstDone || !secondDone;}

    public long next() {
        if (!hasNext()) throw new NoSuchElementException();
        long val = 0;
        if (!firstDone && !secondDone) {
            if (first.peek() <= second.peek()) {
                val = first.next();
                if (!first.hasNext()) {
                    firstDone = true;
                    first.cleanup();
                }
            } else {
                val = second.next();
                if (!second.hasNext()) {
                    secondDone = true;
                    second.cleanup();
                }
            }
        } else if (!firstDone) {
            val = first.next();
            if (!first.hasNext()) {
                firstDone = true;
                first.cleanup();
            }
        } else {
            val = second.next();
            if (!second.hasNext()) {
                secondDone = true;
                second.cleanup();
            }
        }
        return val;
    }

    public void cleanup() {
        if (!firstDone) {
            firstDone = true;
            first.cleanup();
        }
        if (!secondDone) {
            secondDone = true;
            second.cleanup();
        }
    }
}
