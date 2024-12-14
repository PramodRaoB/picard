package picard.sam.markduplicates.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.PeekIterator;

import java.util.*;

/**
 * Utility class that can merge multiple SortingCollections of ReadEndsForMarkDuplicates efficiently
 * using a priority queue.
 */
public class MergingIteratorForSortingCollections<T> implements CloseableIterator<T> {
    private final TreeSet<PeekSortingCollectionIterator> queue;
    List<SortingCollection<T>> collections;
    private final Comparator<T> comparator;
    private boolean isClosed = false;
    private static final Log log = Log.getInstance(SortingCollection.class);

    /**
     * Creates a new MergingSortingCollections instance
     * @param collections List of SortingCollections to merge
     * @param comparator The comparator used by the original SortingCollections
     */
    public MergingIteratorForSortingCollections(List<SortingCollection<T>> collections, Comparator<T> comparator) {
        this.comparator = comparator;
        this.queue = new TreeSet<>(new PeekSortingCollectionIteratorComparator());
        this.collections = collections;
        if (collections.isEmpty()) {
            isClosed = true;
            return;
        }
        initializeQueue(collections);
    }

    private void initializeQueue(List<SortingCollection<T>> collections) {
        int n = 0;
        log.debug(String.format("Creating merging iterator from %d collections", collections.size()));
        for (SortingCollection<T> collection : collections) {
            CloseableIterator<T> iterator = collection.iterator();
            if (iterator.hasNext()) {
                queue.add(new PeekSortingCollectionIterator(iterator, n++));
            } else {
                iterator.close();
            }
        }
    }

    @Override
    public boolean hasNext() {
        return !queue.isEmpty();
    }

    @Override
    public T next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }

        // Get the iterator with the smallest next element
        final PeekSortingCollectionIterator iterator = queue.pollFirst();
        final T result = iterator.next();

        // If this iterator has more elements, add it back to the queue
        if (iterator.hasNext()) {
            queue.add(iterator);
        } else {
            ((CloseableIterator<T>) iterator.getUnderlyingIterator()).close();
        }

        return result;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    @Override
    public void close() {
        if (!isClosed) {
            while (!queue.isEmpty()) {
                PeekSortingCollectionIterator iterator = queue.pollFirst();
                ((CloseableIterator<T>) iterator.getUnderlyingIterator()).close();
            }
            for (SortingCollection<T> collection: this.collections) collection.cleanup();
            isClosed = true;
        }
    }

    /**
     * Wrapper for iterators that allows peeking at the next element
     */
    private class PeekSortingCollectionIterator extends PeekIterator<T> {
        private final int n; // Used for tie-breaking

        PeekSortingCollectionIterator(final Iterator<T> iterator, final int n) {
            super(iterator);
            this.n = n;
        }
    }

    /**
     * Similar to SortingCollection's PeekFileRecordIteratorComparator
     */
    private class PeekSortingCollectionIteratorComparator implements Comparator<PeekSortingCollectionIterator> {
        @Override
        public int compare(PeekSortingCollectionIterator lhs, PeekSortingCollectionIterator rhs) {
            final int result = comparator.compare(lhs.peek(), rhs.peek());
            if (result == 0) return lhs.n - rhs.n;
            else return result;
        }
    }
}