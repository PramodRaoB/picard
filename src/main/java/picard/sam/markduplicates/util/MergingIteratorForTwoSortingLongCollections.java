package picard.sam.markduplicates.util;

import htsjdk.samtools.util.*;

import java.util.*;

/**
 * Very specific utility class that can merge sorting long collections
 * Switch to a priority queue based implementation if required.
 * Currently since its only being used to merge two collections, performance wise it does not make sense
 */
public class MergingIteratorForTwoSortingLongCollections {
    List<SortingLongCollectionForMarkDuplicates> cols;
    List<Boolean> done;
    private static final Log log = Log.getInstance(MergingIteratorForTwoSortingLongCollections.class);

    public MergingIteratorForTwoSortingLongCollections(SortingLongCollectionForMarkDuplicates... args) {
        log.debug(String.format("Creating a merging iterator for Sorting long collections"));
        cols = new ArrayList<>();
        done = new ArrayList<>();
        for (SortingLongCollectionForMarkDuplicates arg: args) {
            if (arg.hasNext()) {
                cols.add(arg);
                done.add(false);
            } else arg.cleanup();
        }
    }

    public boolean hasNext() {return done.stream().anyMatch(n -> n);}

    public long next() {
        if (!hasNext()) throw new NoSuchElementException();
        long val = -1;
        int best = -1;
        for (int idx = 0; idx < cols.size(); idx++) if (!done.get(idx) && (val == -1 || cols.get(idx).peek() < val)) {
            val = cols.get(idx).peek();
            best = idx;
        }
        SortingLongCollectionForMarkDuplicates col = cols.get(best);
        val = col.next();
        if (!col.hasNext()) {
            done.set(best, true);
            col.cleanup();
        }
        return val;
    }

    public void cleanup() {
        for (int idx = 0; idx < cols.size(); idx++) if (!done.get(idx)) {
            done.set(idx, true);
            cols.get(idx).cleanup();
        }
    }
}
