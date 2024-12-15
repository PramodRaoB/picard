/*
 * The MIT License
 *
 * Copyright (c) 2009-2016 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.*;
import picard.sam.util.RepresentativeReadIndexer;

import javax.annotation.Nullable;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.function.Supplier;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = MarkDuplicates.USAGE_SUMMARY + MarkDuplicates.USAGE_DETAILS,
        oneLineSummary = MarkDuplicates.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class MarkDuplicates extends AbstractMarkDuplicatesCommandLineProgram implements MarkDuplicatesHelper {
    static final String USAGE_SUMMARY = "Identifies duplicate reads.  ";
    static final String USAGE_DETAILS = "<p>This tool locates and tags duplicate reads in a SAM, BAM or CRAM file, where duplicate reads are " +
            "defined as originating from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library " +
            "construction using PCR.  See also " +
            "<a href='https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity'>EstimateLibraryComplexity</a>" +
            " for additional notes on PCR duplication artifacts.  Duplicate reads can also result from a single amplification cluster, " +
            "incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.  These duplication artifacts are " +
            "referred to as optical duplicates.</p>" +
            "" +
            "<p>The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file.  " +
            "A BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes.  After duplicate reads are" +
            " collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums " +
            "of their base-quality scores (default method). Note that this is different from directly checking if the sequences match, which " +
            "MarkDuplicates does not do.</p>  " +
            "" +
            "<p>The tool's main output is a new SAM, BAM or CRAM file, in which duplicates have been identified in the SAM flags field for each" +
            " read.  Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.  " +
            "If you are not familiar with this type of annotation, please see the following " +
            "<a href='https://www.broadinstitute.org/gatk/blog?id=7019'>blog post</a> for additional information.</p>" +
            "" +
            "<p>Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of " +
            "duplicate.  To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in  " +
            "the 'optional field' section of a SAM/BAM/CRAM file.  Invoking the TAGGING_POLICY option," +
            " you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no " +
            "duplicates (DontTag).  The records within the output of a SAM/BAM/CRAM file will have values for the 'DT' tag (depending on the invoked " +
            "TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ).  " +
            "This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify " +
            "and differentiate duplicate types.  Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq " +
            "or other data where duplicate sets are extremely large and estimating library complexity is not an aim.  " +
            "Note that without optical duplicate counts, library size estimation will be inaccurate.</p> " +

            "<p>MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.</p>  " +

            "<p>The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different.  " +
            "When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not " +
            "marked as duplicates.  However, when the input is query-sorted (actually query-grouped), " +
            "then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be" +
            " marked as duplicate reads.</p>  " +

            "<p>If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.</p>" +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MarkDuplicates \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=marked_duplicates.bam \\<br />" +
            "      M=marked_dup_metrics.txt" +
            "</pre>" +
            "" +
            "Please see " +
            "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics'>MarkDuplicates</a> " +
            "for detailed explanations of the output metrics." +
            "<hr />";

    /**
     * Enum used to control how duplicates are flagged in the DT optional tag on each read.
     */
    public enum DuplicateTaggingPolicy {
        DontTag, OpticalOnly, All
    }

    /**
     * The optional attribute in SAM/BAM/CRAM files used to store the duplicate type.
     */
    public static final String DUPLICATE_TYPE_TAG = "DT";
    /**
     * The duplicate type tag value for duplicate type: library.
     */
    public static final String DUPLICATE_TYPE_LIBRARY = "LB";
    /**
     * The duplicate type tag value for duplicate type: sequencing (optical & pad-hopping, or "co-localized").
     */
    public static final String DUPLICATE_TYPE_SEQUENCING = "SQ";
    /**
     * The attribute in the SAM/BAM file used to store which read was selected as representative out of a duplicate set
     */
    public static final String DUPLICATE_SET_INDEX_TAG = "DI";
    /**
     * The attribute in the SAM/BAM file used to store the size of a duplicate set
     */
    public static final String DUPLICATE_SET_SIZE_TAG = "DS";

    /**
     * Enum for the possible values that a duplicate read can be tagged with in the DT attribute.
     */
    public enum DuplicateType {
        LIBRARY(DUPLICATE_TYPE_LIBRARY),
        SEQUENCING(DUPLICATE_TYPE_SEQUENCING);

        private final String code;

        DuplicateType(final String code) {
            this.code = code;
        }

        public String code() {
            return this.code;
        }
    }

    private final Log log = Log.getInstance(MarkDuplicates.class);

    /**
     * If more than this many sequences in SAM file, don't spill to disk because there will not
     * be enough file handles.
     */
    @Argument(shortName = "MAX_SEQS",
            doc = "This option is obsolete. ReadEnds will always be spilled to disk.")
    public int MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP = 50000;

    @Argument(shortName = "MAX_FILE_HANDLES",
            doc = "Maximum number of file handles to keep open when spilling read ends to disk. " +
                    "Set this number a little lower than the per-process maximum number of file that may be open. " +
                    "This number can be found by executing the 'ulimit -n' command on a Unix system.")
    public int MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 8000;

    @Argument(doc = "This number, plus the maximum RAM available to the JVM, determine the memory footprint used by " +
            "some of the sorting collections.  If you are running out of memory, try reducing this number.")
    public double SORTING_COLLECTION_SIZE_RATIO = 0.25;

    @Argument(doc = "Barcode SAM tag (ex. BC for 10X Genomics)", optional = true)
    public String BARCODE_TAG = null;

    @Argument(doc = "Read one barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_ONE_BARCODE_TAG = null;

    @Argument(doc = "Read two barcode SAM tag (ex. BX for 10X Genomics)", optional = true)
    public String READ_TWO_BARCODE_TAG = null;

    @Argument(doc = "If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), " +
            "indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two " +
            "reads map to the same portion of the reference only one of which is marked as duplicate. The second " +
            "tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the " +
            "record belongs. This identifier is the index-in-file of the representative read that was selected out " +
            "of the duplicate set.", optional = true)
    public boolean TAG_DUPLICATE_SET_MEMBERS = false;

    @Argument(doc = "If true remove 'optical' duplicates and other duplicates that appear to have arisen from the " +
            "sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. " +
            "If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored.")
    public boolean REMOVE_SEQUENCING_DUPLICATES = false;

    @Argument(doc = "Determines how duplicate types are recorded in the DT optional attribute.")
    public DuplicateTaggingPolicy TAGGING_POLICY = DuplicateTaggingPolicy.DontTag;

    @Argument(doc = "Clear DT tag from input SAM records. Should be set to false if input SAM doesn't have this tag.  Default true")
    public boolean CLEAR_DT = true;

    @Argument(doc = "Treat UMIs as being duplex stranded.  This option requires that the UMI consist of two equal length " +
            "strings that are separated by a hyphen (e.g. 'ATC-GTC'). Reads are considered duplicates if, in addition to standard " +
            "definition, have identical normalized UMIs.  A UMI from the 'bottom' strand is normalized by swapping its content " +
            "around the hyphen (eg. ATC-GTC becomes GTC-ATC).  A UMI from the 'top' strand is already normalized as it is. " +
            "Both reads from a read pair considered top strand if the read 1 unclipped 5' coordinate is less than the read " +
            "2 unclipped 5' coordinate. All chimeric reads and read fragments are treated as having come from the top strand. " +
            "With this option is it required that the BARCODE_TAG hold non-normalized UMIs. Default false.")
    public boolean DUPLEX_UMI = false;

    @Argument(doc = "SAM tag to uniquely identify the molecule from which a read was derived.  Use of this option requires that " +
            "the BARCODE_TAG option be set to a non null value.  Default null.", optional = true)
    public String MOLECULAR_IDENTIFIER_TAG = null;

    @Argument(doc = "Number of reader/writer threads to use. Note that an index file is required for actually multithreading")
    public int NUM_THREADS = 1;

    @Argument(doc = "Optional bedGraph file (e.g. STAR Signal.Unique.str1.out.bg) to use for coverage-aware partitioning", optional = true)
    public String COVERAGE_BEDGRAPH = null;

    @ArgumentCollection
    public MarkDuplicatesForFlowArgumentCollection flowBasedArguments = new MarkDuplicatesForFlowArgumentCollection();

    protected MergingIteratorForSortingCollections<ReadEndsForMarkDuplicates> pairSort;
    protected MergingIteratorForSortingCollections<ReadEndsForMarkDuplicates> fragSort;
    protected GenomicPartitionedCollection<SortingLongCollection> duplicateIndexes;
    protected GenomicPartitionedCollection<SortingLongCollection> opticalDuplicateIndexes;
    protected GenomicPartitionedCollection<SortingCollection<RepresentativeReadIndexer>> representativeReadIndicesForDuplicates;
    protected List<GenomicWindow> windows;
    private List<SortingCollection<ReadEndsForMarkDuplicates>> threadLocalPairSort;
    private List<SortingCollection<ReadEndsForMarkDuplicates>> threadLocalFragSort;
    private List<Set<String>> threadLocalPGIds;
    private Long totalRecords = 0L;
    private Long targetReadsPerChunk;
    private boolean useMultithreading;

    // some calculations are performed using a helper class, which can be parameter specific
    // by default, this instance is the helper
    private MarkDuplicatesHelper calcHelper = this;

    private int numDuplicateIndices = 0;
    static private final long NO_SUCH_INDEX = Long.MAX_VALUE; // needs to be large so that >= test fails for query-sorted traversal

    protected LibraryIdGenerator libraryIdGenerator = null; // this is initialized in buildSortedReadEndLists

    private int getReadOneBarcodeValue(final SAMRecord record) {
        return EstimateLibraryComplexity.getReadBarcodeValue(record, READ_ONE_BARCODE_TAG);
    }

    private int getReadTwoBarcodeValue(final SAMRecord record) {
        return EstimateLibraryComplexity.getReadBarcodeValue(record, READ_TWO_BARCODE_TAG);
    }

    public MarkDuplicates() {
        DUPLICATE_SCORING_STRATEGY = ScoringStrategy.SUM_OF_BASE_QUALITIES;
    }

    public class OperationTimer {
        private static class TimingInfo {
            long startTime;
            long totalTime;
            int calls;

            void start() {
                if (startTime != 0) {
                    throw new IllegalStateException("Timer already started");
                }
                startTime = System.nanoTime();
            }

            void stop() {
                if (startTime == 0) {
                    throw new IllegalStateException("Timer not started");
                }
                totalTime += System.nanoTime() - startTime;
                startTime = 0;
                calls++;
            }

            String getStats() {
                long seconds = totalTime / 1_000_000_000;
                return String.format("%d min %d sec (calls: %d)", seconds / 60, seconds % 60, calls);
            }
        }

        private static final Map<String, TimingInfo> timers = new HashMap<>();

        public static void start(String operation) {
            timers.computeIfAbsent(operation, k -> new TimingInfo()).start();
        }

        public static void stop(String operation) {
            TimingInfo info = timers.get(operation);
            if (info == null) {
                throw new IllegalArgumentException("Unknown operation: " + operation);
            }
            info.stop();
        }

        public static void printStats(Log log) {
            StringBuilder message = new StringBuilder("\nTiming Statistics:\n");
            timers.forEach((op, info) -> message.append(String.format("%s: %s\n", op, info.getStats())));
            log.info(message.toString());
        }
    }

    /**
     * Main work method.  Reads the SAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    protected int doWork() {
        IOUtil.assertInputsAreValid(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        final boolean useBarcodes = (null != BARCODE_TAG || null != READ_ONE_BARCODE_TAG || null != READ_TWO_BARCODE_TAG);

        // use flow based calculation helper?
        if (flowBasedArguments.FLOW_MODE) {
            calcHelper = new MarkDuplicatesForFlowHelper(this);
        }

        reportMemoryStats("Start of doWork");
        log.info("Calculating genomic windows for parallel processing...");
        OperationTimer.start("Calculating genomic windows");
        windows = calculateGenomicWindows();
        OperationTimer.stop("Calculating genomic windows");

        useMultithreading = !windows.isEmpty();
        if (!useMultithreading) {
            log.warn("Could not create processing windows. Proceeding with single-threaded exception.");
        } else {
            log.info("Created " + windows.size() + " processing windows");
        }
        log.info("Reading input file and constructing read end information.");
        OperationTimer.start("Building sorted read end lists");
        buildSortedReadEndLists(useBarcodes);
        if (useMultithreading) windows.sort((w1, w2) -> Long.compare(w2.recordCount, w1.recordCount));
        OperationTimer.stop("Building sorted read end lists");
        reportMemoryStats("After buildSortedReadEndLists");
        OperationTimer.start("Generating duplicate indexes");
        calcHelper.generateDuplicateIndexes(useBarcodes, this.REMOVE_SEQUENCING_DUPLICATES || this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag);
        OperationTimer.stop("Generating duplicate indexes");

        reportMemoryStats("After generateDuplicateIndexes");
        log.info("Marking " + this.numDuplicateIndices + " records as duplicates.");

        if (this.READ_NAME_REGEX == null) {
            log.warn("Skipped optical duplicate cluster discovery; library size estimation may be inaccurate!");
        } else {
            log.info("Found " + (this.libraryIdGenerator.getNumberOfOpticalDuplicateClusters()) + " optical duplicate clusters.");
        }

        for (GenomicWindow window : windows) {
            log.info("Window has " + window.recordCount);
        }

        final SamHeaderAndIterator headerAndIterator = openInputs(false);
        final SAMFileHeader header = headerAndIterator.header;
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();

        final SAMFileHeader outputHeader = header.clone();

        log.info("Reads are assumed to be ordered by: " + sortOrder);

        // Setting the ASSUME_SORT_ORDER to equal queryname is understood to mean that the input is
        // queryname **grouped**. So that's what we set the output order to be, so that the validation will pass
        if (ASSUME_SORT_ORDER == SAMFileHeader.SortOrder.queryname) {
            outputHeader.setGroupOrder(SAMFileHeader.GroupOrder.query);
            outputHeader.setSortOrder(SAMFileHeader.SortOrder.unknown);
            log.info("Output will not be re-sorted. Output header will state SO:unknown GO:query");
        }

        if (ASSUME_SORT_ORDER == null && sortOrder != SAMFileHeader.SortOrder.coordinate && sortOrder != SAMFileHeader.SortOrder.queryname ||
                ASSUME_SORT_ORDER != null && ASSUME_SORT_ORDER != SAMFileHeader.SortOrder.coordinate && ASSUME_SORT_ORDER != SAMFileHeader.SortOrder.queryname) {
            throw new PicardException("This program requires input that are either coordinate or query sorted (according to the header, or at least ASSUME_SORT_ORDER and the content.) " +
                    "Found ASSUME_SORT_ORDER=" + ASSUME_SORT_ORDER + " and header sortorder=" + sortOrder);
        }

        COMMENT.forEach(outputHeader::addComment);

        OperationTimer.start("Main reader/writer loop");
        ExecutorService executor = Executors.newFixedThreadPool(NUM_THREADS);
        List<Future<WindowResult>> futures = new ArrayList<>();

        for (int i = 0; i < windows.size(); i++) {
            final GenomicWindow window = windows.get(i);
            final int windowIndex = i;
            futures.add(executor.submit(() -> processWindow(window, windowIndex)));
        }

        if (!useMultithreading) {
            futures.add(executor.submit(() -> processWindow(null, 0)));
        }

        // Collect and concatenate results
        log.info("Waiting for windows to complete...");
        try {
            List<File> tempFiles = new ArrayList<>(Collections.nCopies(Math.max(windows.size(), 1), null));
            int idx = 0;
            for (Future<WindowResult> future : futures) {
                WindowResult result = future.get();
                final int originalIndex = useMultithreading ? windows.get(idx).windowIndex : 0;
                log.info("Setting index " + originalIndex);
                tempFiles.set(originalIndex, result.tempFile);

                for (Map.Entry<String, DuplicationMetrics> entry :
                        result.metrics.getMetricsByLibraryMap().entrySet()) {
                    DuplicationMetrics globalMetric = this.libraryIdGenerator.getMetricsByLibrary(entry.getKey());
                    if (globalMetric == null) {
                        globalMetric = entry.getValue();
                        this.libraryIdGenerator.addMetricsByLibrary(entry.getKey(), globalMetric);
                    } else {
                        globalMetric.merge(entry.getValue());
                    }
                }
                idx++;
            }

            log.info("Concatenating " + tempFiles.size() + " temp files...");
            concatenateFiles(tempFiles);
        } catch (Exception e) {
            log.error("Exception occurred", e);
        }

        executor.shutdown();
        headerAndIterator.iterator.close();
        OperationTimer.stop("Main reader/writer loop");

        log.info("Closed outputs. Getting more Memory Stats.");
        reportMemoryStats("After output close");

        OperationTimer.printStats(log);

        // Write out the metrics
        finalizeAndWriteMetrics(libraryIdGenerator, getMetricsFile(), METRICS_FILE);

        return 0;
    }

    private WindowResult processWindow(@Nullable GenomicWindow window, int windowIndex) throws IOException {
        long startTime = System.nanoTime();
        File tempFile = File.createTempFile("window_" + windowIndex + "_", ".bam");
        log.info(String.format("Processing starting for window index %d", windowIndex));
        final SamHeaderAndIterator headerAndIterator = openInputs(false);
        final SAMFileHeader header = headerAndIterator.header;
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();

        final SAMFileHeader outputHeader = header.clone();

        // Key: previous PG ID on a SAM Record (or null).  Value: New PG ID to replace it.
        final Map<String, String> chainedPgIds = getChainedPgIds(outputHeader);
        LibraryIdGenerator localMetrics = new LibraryIdGenerator(headerAndIterator.header, flowBasedArguments.FLOW_MODE);

        try (SAMFileWriter out = new SAMFileWriterFactory().makeWriter(outputHeader, true, tempFile, REFERENCE_SEQUENCE)) {
            // Now copy over the file while marking all the necessary indexes as duplicates
            SortingLongCollection duplicateIndexes = this.duplicateIndexes.getPartition(windowIndex);

            SortingLongCollection opticalDuplicateIndexes = this.opticalDuplicateIndexes != null ? this.opticalDuplicateIndexes.getPartition(windowIndex) : null;

            SortingCollection<RepresentativeReadIndexer> representativeReadIndicesForDuplicates = null;
            if (TAG_DUPLICATE_SET_MEMBERS && this.representativeReadIndicesForDuplicates != null) {
                representativeReadIndicesForDuplicates = this.representativeReadIndicesForDuplicates.getPartition(windowIndex);
            }

            long recordInFileIndex = 0;
            long nextOpticalDuplicateIndex = opticalDuplicateIndexes != null && opticalDuplicateIndexes.hasNext() ? opticalDuplicateIndexes.next() : NO_SUCH_INDEX;
            long nextDuplicateIndex = (duplicateIndexes.hasNext() ? duplicateIndexes.next() : NO_SUCH_INDEX);

            // initialize variables for optional representative read tagging
            CloseableIterator<RepresentativeReadIndexer> representativeReadIterator = null;
            RepresentativeReadIndexer rri = null;
            int representativeReadIndexInFile = -1;
            int duplicateSetSize = -1;
            int nextReadInDuplicateSetIndex = -1;
            if (TAG_DUPLICATE_SET_MEMBERS) {
                representativeReadIterator = representativeReadIndicesForDuplicates.iterator();
                if (representativeReadIterator.hasNext()) {
                    rri = representativeReadIterator.next();
                    nextReadInDuplicateSetIndex = rri.readIndexInFile;
                    representativeReadIndexInFile = rri.representativeReadIndexInFile;
                    duplicateSetSize = rri.setSize;
                }
            }

            final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Written");
            CloseableIterator<SAMRecord> iterator;
            if (useMultithreading) {
                final SamReader reader = getReaders(false).get(0);
                iterator = window.isUnmapped() ?
                        reader.queryUnmapped() :
                        reader.query(window.reference, window.start, window.end, false);
                recordInFileIndex = window.startingIndex;
                headerAndIterator.iterator.close();
            }
            else iterator = headerAndIterator.iterator;
            String duplicateQueryName = null;
            String representativeQueryName = null;

            while (iterator.hasNext()) {
                final SAMRecord rec = iterator.next();
                if (window != null && rec.getAlignmentStart() < window.start) continue;

                DuplicationMetrics metrics = AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(rec, header, localMetrics, flowBasedArguments.FLOW_MODE);

                // Now try and figure out the next duplicate index (if going by coordinate. if going by query name, only do this
                // if the query name has changed.
                nextDuplicateIndex = nextIndexIfNeeded(sortOrder, recordInFileIndex, nextDuplicateIndex, duplicateQueryName, rec, duplicateIndexes);

                final boolean isDuplicate = recordInFileIndex == nextDuplicateIndex ||
                        (sortOrder == SAMFileHeader.SortOrder.queryname &&
                                recordInFileIndex > nextDuplicateIndex && rec.getReadName().equals(duplicateQueryName));

                if (isDuplicate) {
                    rec.setDuplicateReadFlag(true);

                    metrics.addDuplicateReadToMetrics(rec);
                } else {
                    rec.setDuplicateReadFlag(false);
                }
                nextOpticalDuplicateIndex = nextIndexIfNeeded(sortOrder, recordInFileIndex, nextOpticalDuplicateIndex, duplicateQueryName, rec, opticalDuplicateIndexes);

                final boolean isOpticalDuplicate = sortOrder == SAMFileHeader.SortOrder.queryname &&
                        recordInFileIndex > nextOpticalDuplicateIndex &&
                        rec.getReadName().equals(duplicateQueryName) ||
                        recordInFileIndex == nextOpticalDuplicateIndex;

                if (CLEAR_DT) {
                    rec.setAttribute(DUPLICATE_TYPE_TAG, null);
                }

                if (this.TAGGING_POLICY != DuplicateTaggingPolicy.DontTag && rec.getDuplicateReadFlag()) {
                    if (isOpticalDuplicate) {
                        rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.SEQUENCING.code());
                    } else if (this.TAGGING_POLICY == DuplicateTaggingPolicy.All) {
                        rec.setAttribute(DUPLICATE_TYPE_TAG, DuplicateType.LIBRARY.code());
                    }
                }

                // Tag any read pair that was in a duplicate set with the duplicate set size and a representative read name
                if (TAG_DUPLICATE_SET_MEMBERS) {
                    final boolean needNextRepresentativeIndex = recordInFileIndex > nextReadInDuplicateSetIndex &&
                            (sortOrder == SAMFileHeader.SortOrder.coordinate || !rec.getReadName().equals(representativeQueryName));
                    if (needNextRepresentativeIndex && representativeReadIterator.hasNext()) {
                        rri = representativeReadIterator.next();
                        nextReadInDuplicateSetIndex = rri.readIndexInFile;
                        representativeReadIndexInFile = rri.representativeReadIndexInFile;
                        duplicateSetSize = rri.setSize;
                    }

                    /* If this record's index is readInDuplicateSetIndex, then it is in a duplicateset.
                    For queryname sorted data, we only have one representativeReadIndex entry per read name, so we need
                    to also look for additional reads with the same name.
                     */
                    final boolean isInDuplicateSet = recordInFileIndex == nextReadInDuplicateSetIndex ||
                            (sortOrder == SAMFileHeader.SortOrder.queryname &&
                                    recordInFileIndex > nextReadInDuplicateSetIndex && rec.getReadName().equals(representativeQueryName));
                    if (isInDuplicateSet) {
                        if (!rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {
                            rec.setAttribute(DUPLICATE_SET_INDEX_TAG, representativeReadIndexInFile);
                            rec.setAttribute(DUPLICATE_SET_SIZE_TAG, duplicateSetSize);
                            representativeQueryName = rec.getReadName();
                        }
                    }
                }

                // Set MOLECULAR_IDENTIFIER_TAG for SAMRecord rec
                if (BARCODE_TAG != null) {
                    UmiUtil.setMolecularIdentifier(rec, "", MOLECULAR_IDENTIFIER_TAG, DUPLEX_UMI);
                }


                // Note, duplicateQueryName must be incremented after we have marked both optical and sequencing duplicates for queryname sorted files.
                if (isDuplicate) {
                    duplicateQueryName = rec.getReadName();
                }

                // Output the record if desired and bump the record index
                recordInFileIndex++;
                if (this.REMOVE_DUPLICATES && rec.getDuplicateReadFlag()) {
                    continue;
                }
                if (this.REMOVE_SEQUENCING_DUPLICATES && isOpticalDuplicate) {
                    continue;
                }
                if (PROGRAM_RECORD_ID != null && pgTagArgumentCollection.ADD_PG_TAG_TO_READS) {
                    rec.setAttribute(SAMTag.PG.name(), chainedPgIds.get(rec.getStringAttribute(SAMTag.PG.name())));
                }
                out.addAlignment(rec);
                progress.record(rec);
            }

            // remember to close the inputs
            iterator.close();

            duplicateIndexes.cleanup();
            if (TAG_DUPLICATE_SET_MEMBERS) {
                representativeReadIndicesForDuplicates.cleanup();
            }
        }
        log.info("For writing window " + windowIndex + " took " + (System.nanoTime() - startTime) / 1_000_000_000 + " seconds.");

        return new WindowResult(tempFile, localMetrics);
    }

    private void concatenateFiles(List<File> tempFiles) throws IOException {
        BamFileIoUtils.gatherWithBlockCopying(tempFiles, OUTPUT, false, false);

        // Clean up temp files
        for (File tempFile : tempFiles) {
            tempFile.delete();
        }
    }


    /**
     * package-visible for testing
     */
    long numOpticalDuplicates() {
        return ((long) this.libraryIdGenerator.getOpticalDuplicatesByLibraryIdMap().getSumOfValues());
    } // cast as long due to returning a double

    /**
     * Print out some quick JVM memory stats.
     */
    private void reportMemoryStats(final String stage) {
        System.gc();
        final Runtime runtime = Runtime.getRuntime();
        log.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
                "; maxMemory: " + runtime.maxMemory());
    }

    /**
     * Calculates genomic windows for parallel processing based on header and index metadata
     */
    private List<GenomicWindow> calculateGenomicWindows() {
        List<GenomicWindow> windows = new ArrayList<>();
        long currentIndex = 0;
        boolean metaNotFound = false;
        int cnt = 0;
        long unmappedCount = 0L;

        try (SamReader reader = SamReaderFactory.makeDefault()
                .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                .open(new File(INPUT.get(0)))) {

            if (reader.hasIndex()) {
                BAMIndex index = reader.indexing().getIndex();
                SAMFileHeader header = reader.getFileHeader();
                for (int i = 0; i < header.getSequenceDictionary().size(); i++) {
                    BAMIndexMetaData metaData = index.getMetaData(i);
                    if (metaData != null) totalRecords += metaData.getAlignedRecordCount();
                    else {
                        log.info("Metadata not found for reference: " + header.getSequence(i).getSequenceName());
                        metaNotFound = true;
                        break;
                    }
                }
                targetReadsPerChunk = totalRecords / NUM_THREADS / 2;

                log.info("Number of reference sequences: " + header.getSequenceDictionary().size());
                CoverageAwareWindowCalculator cov = null;
                if (COVERAGE_BEDGRAPH != null) {
                    // Option to more accurately split windows
                    try {
                        cov = new CoverageAwareWindowCalculator();
                        cov.parseCoverageInformation(COVERAGE_BEDGRAPH);
                    } catch (IOException ignored) {}
                }
                for (int i = 0; i < header.getSequenceDictionary().size(); i++) {
                    SAMSequenceRecord seq = header.getSequence(i);
                    BAMIndexMetaData metadata = index.getMetaData(i);

                    if (metadata != null) {
                        long mappedReads = metadata.getAlignedRecordCount();
                        long unmappedForRef = metadata.getUnalignedRecordCount();
                        if (mappedReads > 0) {
                            GenomicWindow window = new GenomicWindow(
                                    seq.getSequenceName(),
                                    1,
                                    seq.getSequenceLength(),
                                    currentIndex,
                                    mappedReads,
                                    cnt
                            );
                            currentIndex += mappedReads;
                            // Split large windows based on density
                            List<GenomicWindow> splitWindows = splitGenomicWindow(window, cnt, index, i, cov);
                            windows.addAll(splitWindows);
                            cnt += splitWindows.size();
                        }
                        unmappedCount += unmappedForRef;
                    } else {
                        log.info("Metadata not found for reference: " + seq.getSequenceName());
                        metaNotFound = true;
                        break;
                    }
                }

                log.info("Total unmapped count across all references: " + unmappedCount);
                if (unmappedCount > 0) {
                    windows.add(new GenomicWindow(null, -1, -1, currentIndex, unmappedCount, cnt));
                }
            } else metaNotFound = true;
            if (metaNotFound) {
                log.info("Index file is not present or incomplete. Proceeding without multi-threading");
                return new ArrayList<>();
            }
        } catch (IOException ignored) {
            log.info("Index file is not present or incomplete. Proceeding without multi-threading");
            return new ArrayList<>();
        }

        return windows;
    }

    private class CoverageAwareWindowCalculator {
        private class CoverageRegion {
            final int start;
            final int end;
            final double coverage;

            CoverageRegion(int start, int end, double coverage) {
                this.start = start;
                this.end = end;
                this.coverage = coverage;
            }

            long getTotalCoverage() {
                return (long)(coverage * (end - start));
            }
        }

        Map<String, List<CoverageRegion>> coverageByChrom;

        public void parseCoverageInformation(String bedGraphFile) throws IOException {
            this.coverageByChrom = new HashMap<>();

            // Parse bedGraph file
            BufferedReader reader = new BufferedReader(new FileReader(bedGraphFile));
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("track")) continue;

                String[] parts = line.trim().split("\\s+");
                String chrom = parts[0];
                int start = Integer.parseInt(parts[1]);
                int end = Integer.parseInt(parts[2]);
                double coverage = Double.parseDouble(parts[3]);

                if (coverage > 0) {
                    coverageByChrom
                            .computeIfAbsent(chrom, k -> new ArrayList<>())
                            .add(new CoverageRegion(start, end, coverage));
                }
            }
        }

        public List<GenomicWindow> calculateWindows(GenomicWindow window, int windowIndex) {
            List<GenomicWindow> windows = new ArrayList<>();
            long currentIndex = 0;

            List<CoverageRegion> regions = coverageByChrom.getOrDefault(window.reference, Collections.emptyList());
            if (regions.isEmpty()) return windows;

            // Calculate total coverage for this chromosome
            double totalCoverage = regions.stream()
                    .mapToDouble(CoverageRegion::getTotalCoverage)
                    .sum();

            int numWindows = (int) ((double) window.recordCount / targetReadsPerChunk);

            // Calculate target coverage per window
            double targetCoveragePerWindow = totalCoverage / numWindows;

            int genomeStart = regions.get(0).start;

            int currentStart = genomeStart;
            double accumulatedCoverage = 0;
            int regionIndex = 0;

            // Create windows based on coverage
            while (regionIndex < regions.size()) {
                CoverageRegion region = regions.get(regionIndex);
                accumulatedCoverage += region.getTotalCoverage();

                if (accumulatedCoverage >= targetCoveragePerWindow || regionIndex == regions.size() - 1) {
                    // Create window
                    int windowEnd = region.end;
                    windows.add(new MarkDuplicates.GenomicWindow(
                            window.reference,
                            currentStart,
                            windowEnd,
                            currentIndex, // Approximate
                            (long) accumulatedCoverage, // Approximate
                            windowIndex++
                    ));

                    currentStart = windowEnd + 1;
                    currentIndex += accumulatedCoverage;
                    accumulatedCoverage = 0;
                }
                regionIndex++;
            }

            return windows;
        }
    }

    private List<GenomicWindow> splitGenomicWindow(GenomicWindow window, int windowStartIndex, BAMIndex index, int refIndex, @Nullable CoverageAwareWindowCalculator cov) {
        if (cov != null) {
            return cov.calculateWindows(window, windowStartIndex);
        }
        List<GenomicWindow> windows = new ArrayList<>();
        final long TARGET_RECORDS_PER_WINDOW = 3_000_000L;

        try {
            // Get density information
            List<BAMDensityUtil.DensityWindow> densityWindows =
                    BAMDensityUtil.getDensityWindows(index, refIndex, window.end);

            if (densityWindows.isEmpty()) {
                return Collections.singletonList(window);
            }

            // Calculate total density to normalize estimates
            long totalDensity = densityWindows.stream()
                    .mapToLong(dw -> dw.estimatedReads)
                    .sum();

            double readsPerDensityUnit = (double) window.recordCount / totalDensity;

            // Create windows based on density
            long accumulatedReads = 0;
            int currentStart = window.start;
            long currentStartingIndex = window.startingIndex;
            List<BAMDensityUtil.DensityWindow> currentRegions = new ArrayList<>();

            for (BAMDensityUtil.DensityWindow densityWindow : densityWindows) {
                long estimatedReads = (long) (densityWindow.estimatedReads * readsPerDensityUnit);
                accumulatedReads += estimatedReads;
                currentRegions.add(densityWindow);

                if (accumulatedReads >= TARGET_RECORDS_PER_WINDOW) {
                    // Create window from accumulated regions
                    int endPosition = currentRegions.get(currentRegions.size() - 1).end;

                    windows.add(new GenomicWindow(
                            window.reference,
                            currentStart,
                            endPosition,
                            currentStartingIndex,
                            accumulatedReads,
                            windowStartIndex++
                    ));

                    currentStart = endPosition + 1;
                    currentStartingIndex += accumulatedReads;
                    accumulatedReads = 0;
                    currentRegions.clear();
                }
            }

            // Add remaining regions if any
            if (!currentRegions.isEmpty()) {
                int endPosition = currentRegions.get(currentRegions.size() - 1).end;
                windows.add(new GenomicWindow(
                        window.reference,
                        currentStart,
                        endPosition,
                        currentStartingIndex,
                        accumulatedReads,
                        windowStartIndex
                ));
            }

            return windows;

        } catch (Exception e) {
            log.warn("Error splitting window by density, falling back to original window", e);
            return Collections.singletonList(window);
        }
    }

    /**
     * Goes through all the records in a file and generates a set of ReadEndsForMarkDuplicates objects that
     * hold the necessary information (reference sequence, 5' read coordinate) to do
     * duplication, caching to disk as necessary to sort them.
     */
    private void buildSortedReadEndLists(final boolean useBarcodes) {
        final int sizeInBytes;
        if (useBarcodes) {
            sizeInBytes = ReadEndsForMarkDuplicatesWithBarcodes.getSizeOf();
        } else {
            sizeInBytes = ReadEndsForMarkDuplicates.getSizeOf();
        }
        MAX_RECORDS_IN_RAM = (int) (Runtime.getRuntime().maxMemory() / sizeInBytes) / 2;
        int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO) / sizeInBytes);

        final ReadEndsForMarkDuplicatesCodec pairCodec;
        if (useBarcodes) {
            pairCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        } else {
            pairCodec = new ReadEndsForMarkDuplicatesCodec();
        }

        final SamHeaderAndIterator headerAndIterator = openInputs(true);
        final SAMFileHeader header = headerAndIterator.header;
        ReadEndsForMarkDuplicatesCodec diskCodec;
        if (useBarcodes) diskCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        else diskCodec = new ReadEndsForMarkDuplicatesCodec();
        final ReadEndsForMarkDuplicatesMap tmp = useMultithreading ? new DiskBasedReadEndsForMarkDuplicatesMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, diskCodec) : null;
        this.threadLocalFragSort = new ArrayList<>(Collections.nCopies(NUM_THREADS, null));
        this.threadLocalPairSort = new ArrayList<>(Collections.nCopies(NUM_THREADS, null));
        this.threadLocalPGIds = new ArrayList<>(Collections.nCopies(NUM_THREADS, null));
        List<SingleMemoryBasedReadEndsForMarkDuplicatesMap> extraList = new ArrayList<>();

        if (null == this.libraryIdGenerator) {
            this.libraryIdGenerator = new LibraryIdGenerator(header, flowBasedArguments.FLOW_MODE);
        }

        ExecutorService executor = Executors.newFixedThreadPool(NUM_THREADS);
        List<Future<ReadEndsForMarkDuplicatesMap>> futures = new ArrayList<>();

        for (GenomicWindow window : windows) {
            futures.add(executor.submit(() -> processWindowForReadEnds(window, useBarcodes)));
        }

        if (!useMultithreading) futures.add(executor.submit(() -> processWindowForReadEnds(null, useBarcodes)));
        headerAndIterator.iterator.close();

        SortingCollection<ReadEndsForMarkDuplicates> pairSort = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                pairCodec,
                new ReadEndsMDComparator(useBarcodes),
                maxInMemory,
                TMP_DIR);
        try {
            for (Future<ReadEndsForMarkDuplicatesMap> future : futures) {
                if (useMultithreading) extraList.add((SingleMemoryBasedReadEndsForMarkDuplicatesMap) future.get());
            }
        } catch (Exception e) {
            throw new PicardException("Error processing windows", e);
        }
        if (useMultithreading) {
            long startingIndex = 0;
            long processedSequentially = 0;
            long cnt = 0;
            for (int idx = 0; idx < windows.size(); idx++) cnt += extraList.get(idx).size();
            log.info("Going to process " + cnt + " records sequentially");
            for (int idx = 0; idx < windows.size(); idx++) {
                log.info("Starting index " + idx + " with size: " + extraList.get(idx).size());
                // Set the correct starting index for later steps
                windows.get(idx).setStartingIndex(startingIndex);
                SingleMemoryBasedReadEndsForMarkDuplicatesMap map = extraList.get(idx);
                for (Map.Entry<String, ReadEndsForMarkDuplicates> entry : map.entrySet()) {
                    processedSequentially++;
                    String key = entry.getKey();
                    ReadEndsForMarkDuplicates fragmentEnd = entry.getValue();
                    // Correct the index of this readEnd
                    fragmentEnd.read1IndexInFile = fragmentEnd.read1IndexInFile + startingIndex;
                    ReadEndsForMarkDuplicates pairedEnds = tmp.remove(fragmentEnd.read1ReferenceIndex, key);
                    // See if we've already seen the first end or not
                    if (pairedEnds == null) {
                        // at this point pairedEnds and fragmentEnd are the same, but we need to make
                        // a copy since pairedEnds will be modified when the mate comes along.
                        pairedEnds = fragmentEnd.clone();
                        tmp.put(pairedEnds.read2ReferenceIndex, key.toString(), pairedEnds);
                    } else {
                        final int matesRefIndex = fragmentEnd.read1ReferenceIndex;
                        final int matesCoordinate = fragmentEnd.read1Coordinate;

                        // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
                        // before updating the orientation later.
                        if (fragmentEnd.firstOfFlag == 1) {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(fragmentEnd.orientation == ReadEnds.R, pairedEnds.orientation == ReadEnds.R);
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readOneBarcode = ((ReadEndsForMarkDuplicatesWithBarcodes) fragmentEnd).readOneBarcode;
                            }
                        } else {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, fragmentEnd.orientation == ReadEnds.R);
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readOneBarcode = ((ReadEndsForMarkDuplicatesWithBarcodes) fragmentEnd).readTwoBarcode;
                            }
                        }

                        // If the other read is actually later, simply add the other read's data as read2, else flip the reads
                        if (matesRefIndex > pairedEnds.read1ReferenceIndex ||
                                (matesRefIndex == pairedEnds.read1ReferenceIndex && matesCoordinate >= pairedEnds.read1Coordinate)) {
                            pairedEnds.read2ReferenceIndex = matesRefIndex;
                            pairedEnds.read2Coordinate = matesCoordinate;
                            pairedEnds.read2IndexInFile = fragmentEnd.read1IndexInFile;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, fragmentEnd.orientation == ReadEnds.R);

                            // if the two read ends are in the same position, pointing in opposite directions,
                            // the orientation is undefined and the procedure above
                            // will depend on the order of the reads in the file.
                            // To avoid this, we set it explicitly (to FR):
                            if (pairedEnds.read2ReferenceIndex == pairedEnds.read1ReferenceIndex &&
                                    pairedEnds.read2Coordinate == pairedEnds.read1Coordinate &&
                                    pairedEnds.orientation == ReadEnds.RF) {
                                pairedEnds.orientation = ReadEnds.FR;
                            }
                        } else {
                            pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
                            pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
                            pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                            pairedEnds.read1ReferenceIndex = matesRefIndex;
                            pairedEnds.read1Coordinate = matesCoordinate;
                            pairedEnds.read1IndexInFile = fragmentEnd.read1IndexInFile;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(fragmentEnd.orientation == ReadEnds.R,
                                    pairedEnds.orientation == ReadEnds.R);
                        }

                        // pairedEnds.score += calcHelper.getReadDuplicateScore(rec, pairedEnds);
                        pairedEnds.score += fragmentEnd.score;
                        pairSort.add(pairedEnds);
                    }
                }
                startingIndex += windows.get(idx).recordCount;
            }
            log.info("Processed " + processedSequentially + " records sequentially.");
        }
        this.threadLocalPairSort.add(pairSort);
        if (useMultithreading) log.info(tmp.size() + " pairs never matched.");
        log.info("Merging built read end lists");
        for (SortingCollection<ReadEndsForMarkDuplicates> col : this.threadLocalPairSort) if (col != null) col.doneAdding();
        for (SortingCollection<ReadEndsForMarkDuplicates> col : this.threadLocalFragSort) if (col != null) col.doneAdding();
        for (Set<String> st : this.threadLocalPGIds) if (st != null) mergeCollection(this.pgIdsSeen, st);
        this.fragSort = new MergingIteratorForSortingCollections<>(threadLocalFragSort, new ReadEndsMDComparator(useBarcodes));
        this.pairSort = new MergingIteratorForSortingCollections<>(threadLocalPairSort, new ReadEndsMDComparator(useBarcodes));

        executor.shutdown();
    }

    private ReadEndsForMarkDuplicatesMap processWindowForReadEnds(@Nullable GenomicWindow window, boolean useBarcodes) {
        long startTime = System.nanoTime();
        final int sizeInBytes;
        if (useBarcodes) {
            sizeInBytes = ReadEndsForMarkDuplicatesWithBarcodes.getSizeOf();
        } else {
            sizeInBytes = ReadEndsForMarkDuplicates.getSizeOf();
        }
        int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * SORTING_COLLECTION_SIZE_RATIO) / sizeInBytes);
        if (useMultithreading) maxInMemory /= windows.size();
        log.info("Will retain up to " + maxInMemory + " data points before spilling to disk.");

        final ReadEndsForMarkDuplicatesCodec fragCodec, pairCodec, diskCodec;
        if (useBarcodes) {
            fragCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
            pairCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
            diskCodec = new ReadEndsForMarkDuplicatesWithBarcodesCodec();
        } else {
            fragCodec = new ReadEndsForMarkDuplicatesCodec();
            pairCodec = new ReadEndsForMarkDuplicatesCodec();
            diskCodec = new ReadEndsForMarkDuplicatesCodec();
        }

        int tid = (int) (Thread.currentThread().getId() % NUM_THREADS);

        SortingCollection<ReadEndsForMarkDuplicates> pairSort = this.threadLocalPairSort.get(tid);
        if (pairSort == null) pairSort = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                pairCodec,
                new ReadEndsMDComparator(useBarcodes),
                maxInMemory,
                TMP_DIR);

        SortingCollection<ReadEndsForMarkDuplicates> fragSort = this.threadLocalPairSort.get(tid);
        if (fragSort == null) fragSort = SortingCollection.newInstance(ReadEndsForMarkDuplicates.class,
                fragCodec,
                new ReadEndsMDComparator(useBarcodes),
                maxInMemory,
                TMP_DIR);

        Set<String> pgIdsSeen = this.threadLocalPGIds.get(tid);
        if (pgIdsSeen == null) pgIdsSeen = new HashSet<>();

        final SamHeaderAndIterator headerAndIterator = openInputs(true);
        final SAMFileHeader.SortOrder assumedSortOrder = headerAndIterator.header.getSortOrder();
        final SAMFileHeader header = headerAndIterator.header;
        final ReadEndsForMarkDuplicatesMap tmp =
                (useMultithreading || assumedSortOrder == SAMFileHeader.SortOrder.queryname)
                        ? new SingleMemoryBasedReadEndsForMarkDuplicatesMap()
                        : new DiskBasedReadEndsForMarkDuplicatesMap(MAX_FILE_HANDLES_FOR_READ_ENDS_MAP, diskCodec);

        long index = 0;
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e6, "Read");
        final CloseableIterator<SAMRecord> iterator;
        if (useMultithreading) {
            final SamReader reader = getReaders(false).get(0);
            iterator = window.isUnmapped() ?
                    reader.queryUnmapped() :
                    reader.query(window.reference, window.start, window.end, false);
            headerAndIterator.iterator.close();
        }
        else iterator = headerAndIterator.iterator;

        String duplicateQueryName = null;
        long duplicateIndex = NO_SUCH_INDEX;

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            if (window != null && rec.getAlignmentStart() < window.start) continue;

            // This doesn't have anything to do with building sorted ReadEnd lists, but it can be done in the same pass
            // over the input
            if (PROGRAM_RECORD_ID != null) {
                // Gather all PG IDs seen in merged input files in first pass.  These are gathered for two reasons:
                // - to know how many different PG records to create to represent this program invocation.
                // - to know what PG IDs are already used to avoid collisions when creating new ones.
                // Note that if there are one or more records that do not have a PG tag, then a null value
                // will be stored in this set.
                pgIdsSeen.add(rec.getStringAttribute(SAMTag.PG.name()));
            }

            // If working in query-sorted, need to keep index of first record with any given query-name.
            if (assumedSortOrder == SAMFileHeader.SortOrder.queryname && !rec.getReadName().equals(duplicateQueryName)) {
                duplicateQueryName = rec.getReadName();
                duplicateIndex = index;
            }

            if (rec.getReadUnmappedFlag()) {
                if (rec.getReferenceIndex() == -1 && assumedSortOrder == SAMFileHeader.SortOrder.coordinate) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue (only in coordinate sort).
                    break;
                }
                // If this read is unmapped but sorted with the mapped reads, just skip it.

            } else if (!rec.isSecondaryOrSupplementary()) {
                final long indexForRead = assumedSortOrder == SAMFileHeader.SortOrder.queryname ? duplicateIndex : index;
                final ReadEndsForMarkDuplicates fragmentEnd = calcHelper.buildReadEnds(header, indexForRead, rec, useBarcodes);
                fragSort.add(fragmentEnd);

                if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                    final StringBuilder key = new StringBuilder();
                    key.append(rec.getReadGroup().getReadGroupId());
                    key.append(rec.getReadName());
                    ReadEndsForMarkDuplicates pairedEnds = tmp.remove(rec.getReferenceIndex(), key.toString());

                    // See if we've already seen the first end or not
                    if (pairedEnds == null) {
                        // at this point pairedEnds and fragmentEnd are the same, but we need to make
                        // a copy since pairedEnds will be modified when the mate comes along.
                        pairedEnds = fragmentEnd.clone();
                        if (useBarcodes) {
                            ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readOneBarcode = getReadOneBarcodeValue(rec);
                            ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readTwoBarcode = getReadTwoBarcodeValue(rec);
                        }
                        tmp.put(pairedEnds.read2ReferenceIndex, key.toString(), pairedEnds);
                    } else {
                        final int matesRefIndex = fragmentEnd.read1ReferenceIndex;
                        final int matesCoordinate = fragmentEnd.read1Coordinate;

                        // Set orientationForOpticalDuplicates, which always goes by the first then the second end for the strands.  NB: must do this
                        // before updating the orientation later.
                        if (rec.getFirstOfPairFlag()) {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(), pairedEnds.orientation == ReadEnds.R);
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readOneBarcode = getReadOneBarcodeValue(rec);
                            }
                        } else {
                            pairedEnds.orientationForOpticalDuplicates = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R, rec.getReadNegativeStrandFlag());
                            if (useBarcodes) {
                                ((ReadEndsForMarkDuplicatesWithBarcodes) pairedEnds).readTwoBarcode = getReadTwoBarcodeValue(rec);
                            }
                        }

                        // If the other read is actually later, simply add the other read's data as read2, else flip the reads
                        if (matesRefIndex > pairedEnds.read1ReferenceIndex ||
                                (matesRefIndex == pairedEnds.read1ReferenceIndex && matesCoordinate >= pairedEnds.read1Coordinate)) {
                            pairedEnds.read2ReferenceIndex = matesRefIndex;
                            pairedEnds.read2Coordinate = matesCoordinate;
                            pairedEnds.read2IndexInFile = indexForRead;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(pairedEnds.orientation == ReadEnds.R,
                                    rec.getReadNegativeStrandFlag());

                            // if the two read ends are in the same position, pointing in opposite directions,
                            // the orientation is undefined and the procedure above
                            // will depend on the order of the reads in the file.
                            // To avoid this, we set it explicitly (to FR):
                            if (pairedEnds.read2ReferenceIndex == pairedEnds.read1ReferenceIndex &&
                                    pairedEnds.read2Coordinate == pairedEnds.read1Coordinate &&
                                    pairedEnds.orientation == ReadEnds.RF) {
                                pairedEnds.orientation = ReadEnds.FR;
                            }
                        } else {
                            pairedEnds.read2ReferenceIndex = pairedEnds.read1ReferenceIndex;
                            pairedEnds.read2Coordinate = pairedEnds.read1Coordinate;
                            pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                            pairedEnds.read1ReferenceIndex = matesRefIndex;
                            pairedEnds.read1Coordinate = matesCoordinate;
                            pairedEnds.read1IndexInFile = indexForRead;
                            pairedEnds.orientation = ReadEnds.getOrientationByte(rec.getReadNegativeStrandFlag(),
                                    pairedEnds.orientation == ReadEnds.R);
                        }

                        pairedEnds.score += fragmentEnd.score;
                        pairSort.add(pairedEnds);
                    }
                }
            }

            // Print out some stats every 1m reads
            ++index;
            if (progress.record(rec)) {
                log.info("Tracking " + tmp.size() + " as yet unmatched pairs. " + tmp.sizeInRam() + " records in RAM.");
            }
        }

        // Set the correct record count for the current partition
        if (window != null) window.setRecordCount(index);

        log.info("Read " + index + " records. " + tmp.size() + " pairs never matched.");
        iterator.close();

        log.info("For building " + index + " records - took " + (System.nanoTime() - startTime) / 1_000_000_000 + " seconds.");
        return tmp;
    }

    private <T> void mergeCollection(SortingCollection<T> out, SortingCollection<T> in) {
        for (final T v: in) out.add(v);
    }
    private <T> void mergeCollection(Set<T> out, Set<T> in) {
        out.addAll(in);
    }

    /**
     * Calculates score for the duplicate read
     * @param rec - read
     * @param pairedEnds - location of the read ends
     * @return - read score calculated according to the DUPLICATE_SCORING_STRATEGY:
     *         SUM_OF_BASE_QUALITIES, (default)
     *         TOTAL_MAPPED_REFERENCE_LENGTH,
     *         RANDOM
     *
     */
    public short getReadDuplicateScore(final SAMRecord rec, final ReadEndsForMarkDuplicates pairedEnds) {
        return DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
    }

    /**
     * Builds a read ends object that represents a single read.
     */
    public ReadEndsForMarkDuplicates buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec, final boolean useBarcodes) {
        final ReadEndsForMarkDuplicates ends;

        if (useBarcodes) {
            ends = new ReadEndsForMarkDuplicatesWithBarcodes();
        } else {
            ends = new ReadEndsForMarkDuplicates();
        }
        ends.read1ReferenceIndex = rec.getReferenceIndex();
        ends.read1Coordinate = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.read1IndexInFile = index;
        ends.score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);

        if (rec.getFirstOfPairFlag()) ends.firstOfFlag = 1;

        // Doing this lets the ends object know that it's part of a pair
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2ReferenceIndex = rec.getMateReferenceIndex();
        }

        // Fill in the library ID
        ends.libraryId = libraryIdGenerator.getLibraryId(rec);

        // Fill in the location information for optical duplicates
        if (this.opticalDuplicateFinder.addLocationInformation(rec.getReadName(), ends)) {
            // calculate the RG number (nth in list)
            ends.readGroup = 0;
            final String rg = (String) rec.getAttribute(ReservedTagConstants.READ_GROUP_ID);
            final List<SAMReadGroupRecord> readGroups = header.getReadGroups();

            if (rg != null && readGroups != null) {
                for (final SAMReadGroupRecord readGroup : readGroups) {
                    if (readGroup.getReadGroupId().equals(rg)) {
                        break;
                    } else {
                        ends.readGroup++;
                    }
                }
            }
        }

        if (useBarcodes) {
            final ReadEndsForMarkDuplicatesWithBarcodes endsWithBarcode = (ReadEndsForMarkDuplicatesWithBarcodes) ends;
            final String topStrandNormalizedUmi = UmiUtil.getTopStrandNormalizedUmi(rec, BARCODE_TAG, DUPLEX_UMI);
            endsWithBarcode.barcode = Objects.hash(topStrandNormalizedUmi);

            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                endsWithBarcode.readOneBarcode = getReadOneBarcodeValue(rec);
            } else {
                endsWithBarcode.readTwoBarcode = getReadTwoBarcodeValue(rec);
            }
        }

        return ends;
    }

    /**
     * Goes through the accumulated ReadEndsForMarkDuplicates objects and determines which of them are
     * to be marked as duplicates.
     */
    protected void sortIndicesForDuplicates(final boolean indexOpticalDuplicates){
        final int entryOverhead;
        if (TAG_DUPLICATE_SET_MEMBERS) {
            // Memory requirements for RepresentativeReadIndexer:
            // three int entries + overhead: (3 * 4) + 4 = 16 bytes
            entryOverhead = 16;
        } else {
            entryOverhead = SortingLongCollection.SIZEOF;
        }
        // Keep this number from getting too large even if there is a huge heap.
        int maxInMemory = (int) Math.min((Runtime.getRuntime().maxMemory() * 0.25) / entryOverhead, (double) (Integer.MAX_VALUE - 5));
        if (useMultithreading) maxInMemory /= windows.size();
        // If we're also tracking optical duplicates, reduce maxInMemory, since we'll need two sorting collections
        if (indexOpticalDuplicates) {
            maxInMemory /= ((entryOverhead + SortingLongCollection.SIZEOF) / entryOverhead);
            int finalMaxInMemory = maxInMemory;
            this.opticalDuplicateIndexes = new GenomicPartitionedCollection<>(windows, () -> new SortingLongCollection(finalMaxInMemory, TMP_DIR.toArray(new File[0])));
        }

        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        int finalMaxInMemory = maxInMemory;
        this.duplicateIndexes = new GenomicPartitionedCollection<>(windows, () -> new SortingLongCollection(finalMaxInMemory, TMP_DIR.toArray(new File[0])));
        if (TAG_DUPLICATE_SET_MEMBERS) {
            final RepresentativeReadIndexerCodec representativeIndexCodec = new RepresentativeReadIndexerCodec();
            this.representativeReadIndicesForDuplicates =
                    new GenomicPartitionedCollection<>(windows,
                            () -> SortingCollection.newInstance(
                                    RepresentativeReadIndexer.class,
                                    representativeIndexCodec,
                                    Comparator.comparing(read -> read.readIndexInFile),
                                    finalMaxInMemory,
                                    TMP_DIR
                            ));
        }
    }

    public void generateDuplicateIndexes(final boolean useBarcodes, final boolean indexOpticalDuplicates) {
        sortIndicesForDuplicates(indexOpticalDuplicates);

        ReadEndsForMarkDuplicates firstOfNextChunk = null;
        final List<ReadEndsForMarkDuplicates> nextChunk = new ArrayList<>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        while (this.pairSort.hasNext()) {
            final ReadEndsForMarkDuplicates next = this.pairSort.next();
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, true, useBarcodes)) {
                nextChunk.add(next);
            } else {
                handleChunk(nextChunk);
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
            }
        }
        handleChunk(nextChunk);

        this.pairSort.close();
        this.pairSort = null;

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        firstOfNextChunk = null;

        while (this.fragSort.hasNext()) {
            final ReadEndsForMarkDuplicates next = this.fragSort.next();
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false, useBarcodes)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
            } else {
                if (nextChunk.size() > 1 && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }
                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);

        this.fragSort.close();
        this.fragSort = null;

        log.info("Sorting list of duplicate records.");
        this.duplicateIndexes.doneAdding();
        if (this.opticalDuplicateIndexes != null) {
            this.opticalDuplicateIndexes.doneAdding();
        }
        if (TAG_DUPLICATE_SET_MEMBERS) {
            this.representativeReadIndicesForDuplicates.doneAdding();
        }
    }

    protected void handleChunk(List<ReadEndsForMarkDuplicates> nextChunk) {
        if (nextChunk.size() > 1) {
            markDuplicatePairs(nextChunk);
            if (TAG_DUPLICATE_SET_MEMBERS) {
                addRepresentativeReadIndex(nextChunk);
            }
        } else if (nextChunk.size() == 1) {
            addSingletonToCount(libraryIdGenerator);
        }
    }

    protected boolean areComparableForDuplicates(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs, final boolean compareRead2, final boolean useBarcodes) {
        boolean areComparable = lhs.libraryId == rhs.libraryId;

        if (useBarcodes && areComparable) { // areComparable is useful here to avoid the casts below
            final ReadEndsForMarkDuplicatesWithBarcodes lhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) lhs;
            final ReadEndsForMarkDuplicatesWithBarcodes rhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) rhs;
            areComparable = lhsWithBarcodes.barcode == rhsWithBarcodes.barcode &&
                    lhsWithBarcodes.readOneBarcode == rhsWithBarcodes.readOneBarcode &&
                    lhsWithBarcodes.readTwoBarcode == rhsWithBarcodes.readTwoBarcode;
        }

        if (areComparable) {
            areComparable = lhs.read1ReferenceIndex == rhs.read1ReferenceIndex &&
                    lhs.read1Coordinate == rhs.read1Coordinate &&
                    lhs.orientation == rhs.orientation;
        }

        if (areComparable && compareRead2) {
            areComparable = lhs.read2ReferenceIndex == rhs.read2ReferenceIndex &&
                    lhs.read2Coordinate == rhs.read2Coordinate;
        }

        return areComparable;
    }

    private void addIndexAsDuplicate(final long bamIndex) {
        this.duplicateIndexes.add(bamIndex);
        ++this.numDuplicateIndices;
    }

    private void addRepresentativeReadOfDuplicateSet(final long representativeReadIndexInFile, final int setSize, final long read1IndexInFile) {
        final RepresentativeReadIndexer rri = new RepresentativeReadIndexer();
        rri.representativeReadIndexInFile = (int) representativeReadIndexInFile;
        rri.setSize = setSize;
        rri.readIndexInFile = (int) read1IndexInFile;
        this.representativeReadIndicesForDuplicates.add(rri);
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and identify the representative read based on
     * quality score. For all members of the duplicate set, add the read1 index-in-file of the representative
     * read to the records of the first and second in a pair. This value becomes is used for
     * the 'DI' tag.
     */
    private void addRepresentativeReadIndex(final List<ReadEndsForMarkDuplicates> list) {
        short maxScore = 0;
        ReadEndsForMarkDuplicates best = null;

        /** All read ends should have orientation FF, FR, RF, or RR **/
        for (final ReadEndsForMarkDuplicates end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        // for read name (for representative read name), add the last of the pair that was examined
        for (final ReadEndsForMarkDuplicates end : list) {
            addRepresentativeReadOfDuplicateSet(best.read1IndexInFile, list.size(), end.read1IndexInFile);
            if (end.read1IndexInFile != end.read2IndexInFile) {
                addRepresentativeReadOfDuplicateSet(best.read1IndexInFile, list.size(), end.read2IndexInFile);
            }
        }
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This assumes that the list contains objects representing pairs.
     */
    private void markDuplicatePairs(final List<ReadEndsForMarkDuplicates> list) {
        short maxScore = 0;
        ReadEndsForMarkDuplicates best = null;

        /** All read ends should have orientation FF, FR, RF, or RR **/
        for (final ReadEndsForMarkDuplicates end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        if (this.READ_NAME_REGEX != null) {
            AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(list, best, opticalDuplicateFinder, libraryIdGenerator);
        }

        for (final ReadEndsForMarkDuplicates end : list) {
            if (end != best) {
                addIndexAsDuplicate(end.read1IndexInFile);

                // in query-sorted case, these will be the same.
                // TODO: also in coordinate sorted, when one read is unmapped
                if (end.read2IndexInFile != end.read1IndexInFile) {
                    addIndexAsDuplicate(end.read2IndexInFile);
                }

                if (end.isOpticalDuplicate && this.opticalDuplicateIndexes != null) {
                    this.opticalDuplicateIndexes.add(end.read1IndexInFile);
                    // We expect end.read2IndexInFile==read1IndexInFile when we are in queryname sorted files, as the read-pairs
                    // will be sorted together and nextIndexIfNeeded() will only pull one index from opticalDuplicateIndexes.
                    // This means that in queryname sorted order we will only pull from the sorting collection once,
                    // where as we would pull twice for coordinate sorted files.
                    if (end.read2IndexInFile != end.read1IndexInFile) {
                        this.opticalDuplicateIndexes.add(end.read2IndexInFile);
                    }
                }
            }
        }
    }

    /**
     * Method for deciding when to pull from the SortingLongCollection for the next read based on sorting order.
     * - If file is queryname sorted then we expect one index per pair of reads, so we only want to iterate when we
     * are no longer reading from that read-pair.
     * - If file is coordinate-sorted we want to base our iteration entirely on the indexes of both reads in the pair
     * <p>
     * This logic is applied to both Optical and Library duplicates
     *
     * @param sortOrder          Sort order for the underlying bam file
     * @param recordInFileIndex  Index of the current sam record rec
     * @param nextDuplicateIndex Index of next expected duplicate (optical or otherwise) in the file
     * @param lastQueryName      Name of the last read seen (for keeping queryname sorted groups together)
     * @param rec                Current record to compare against
     * @param duplicateIndexes   DuplicateIndexes collection to iterate over
     * @return the duplicate after iteration
     */
    private long nextIndexIfNeeded(final SAMFileHeader.SortOrder sortOrder, final long recordInFileIndex, final long nextDuplicateIndex, final String lastQueryName, final SAMRecord rec, final SortingLongCollection duplicateIndexes) {
        // Manage the flagging of optical/sequencing duplicates
        // Possibly figure out the next opticalDuplicate index (if going by coordinate, if going by query name, only do this
        // if the query name has changed)
        final boolean needNextDuplicateIndex = recordInFileIndex > nextDuplicateIndex &&
                (sortOrder == SAMFileHeader.SortOrder.coordinate || !rec.getReadName().equals(lastQueryName));

        if (needNextDuplicateIndex) {
            return (duplicateIndexes.hasNext() ? duplicateIndexes.next() : NO_SUCH_INDEX);
        }
        return nextDuplicateIndex;
    }

    /**
     * Takes a list of ReadEndsForMarkDuplicates objects and removes from it all objects that should
     * not be marked as duplicates.  This will set the duplicate index for only list items are fragments.
     *
     * @param containsPairs true if the list also contains objects containing pairs, false otherwise.
     */
    protected void markDuplicateFragments(final List<ReadEndsForMarkDuplicates> list, final boolean containsPairs) {
        if (containsPairs) {
            for (final ReadEndsForMarkDuplicates end : list) {
                if (!end.isPaired()) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        } else {
            short maxScore = 0;
            ReadEndsForMarkDuplicates best = null;
            for (final ReadEndsForMarkDuplicates end : list) {
                if (end.score > maxScore || best == null) {
                    maxScore = end.score;
                    best = end;
                }
            }

            for (final ReadEndsForMarkDuplicates end : list) {
                if (end != best) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        }
    }

    static class GenomicWindow {
        final String reference;
        final int start;
        final int end;
        long startingIndex;
        long recordCount;
        int windowIndex;

        GenomicWindow(String reference, int start, int end, long startingIndex, long recordCount, int windowIndex) {
            this.reference = reference;
            this.start = start;
            this.end = end;
            this.startingIndex = startingIndex;
            this.recordCount = recordCount;
            this.windowIndex = windowIndex;
        }

        void setRecordCount(long recordCount) {
            this.recordCount = recordCount;
        }

        boolean isUnmapped() {
            return reference == null;
        }

        public void setStartingIndex(long startingIndex) {
            this.startingIndex = startingIndex;
        }
    }

    private class WindowResult {
        final File tempFile;
        final LibraryIdGenerator metrics;

        WindowResult(File tempFile, LibraryIdGenerator metrics) {
            this.tempFile = tempFile;
            this.metrics = metrics;
        }
    }

    private class ThreadLocalReadEnds {
        final SortingCollection<ReadEndsForMarkDuplicates> localPairSort;
        final SortingCollection<ReadEndsForMarkDuplicates> localFragSort;
        final ReadEndsForMarkDuplicatesMap map;
        Set<String> pgIdsSeen;

        ThreadLocalReadEnds(SortingCollection<ReadEndsForMarkDuplicates> pairSort,
                            SortingCollection<ReadEndsForMarkDuplicates> fragSort,
                            ReadEndsForMarkDuplicatesMap map,
                            Set<String> pgIdsSeen) {
            this.localPairSort = pairSort;
            this.localFragSort = fragSort;
            this.map = map;
            this.pgIdsSeen = pgIdsSeen;
        }
    }

    class GenomicPartitionedCollection<T> {
        private final List<T> collections;
        private final List<GenomicWindow> windows;

        public GenomicPartitionedCollection(List<GenomicWindow> windows, Supplier<T> collectionFactory) {
            this.windows = windows;
            this.collections = new ArrayList<>(Math.max(1, windows.size()));

            // Create one collection per window, or just one if no windows
            for (int i = 0; i < Math.max(1, windows.size()); i++) {
                collections.add(collectionFactory.get());
            }
        }

        public void add(Object value) {
            long fileIndex;
            if (value instanceof Long) {
                fileIndex = (Long) value;
            } else if (value instanceof RepresentativeReadIndexer) {
                fileIndex = ((RepresentativeReadIndexer) value).readIndexInFile;
            } else {
                throw new IllegalArgumentException("Unsupported value type");
            }

            int partition = findPartition(fileIndex);
            T collection = collections.get(partition);

            // Handle different collection types
            if (collection instanceof SortingLongCollection) {
                ((SortingLongCollection) collection).add((Long) value);
            } else if (collection instanceof SortingCollection) {
                ((SortingCollection) collection).add(value);
            }
        }

        public T getPartition(int index) {
            return collections.get(index);
        }

        public void doneAdding() {
            for (T collection : collections) {
                if (collection instanceof SortingLongCollection) {
                    ((SortingLongCollection) collection).doneAddingStartIteration();
                } else if (collection instanceof SortingCollection) {
                    ((SortingCollection) collection).doneAdding();
                }
            }
        }

        private int findPartition(long fileIndex) {
            if (windows.isEmpty()) return 0;

            for (int i = 0; i < windows.size(); i++) {
                GenomicWindow window = windows.get(i);
                if (fileIndex >= window.startingIndex &&
                        fileIndex < window.startingIndex + window.recordCount) {
                    return i;
                }
            }
            throw new IllegalArgumentException("File index " + fileIndex +
                    " does not belong to any window");
        }
    }

    class BAMDensityUtil {
        // Represents density information for a genomic region
        static class DensityWindow {
            final int start;
            final int end;
            final long estimatedReads;

            DensityWindow(int start, int end, long estimatedReads) {
                this.start = start;
                this.end = end;
                this.estimatedReads = estimatedReads;
            }
        }

        // Get density windows for a reference sequence using level 2 bins (1Mb regions)
        static List<DensityWindow> getDensityWindows(BAMIndex index, int referenceIndex, int sequenceLength) {
            List<DensityWindow> densityWindows = new ArrayList<>();

            if (!(index instanceof BrowseableBAMIndex)) {
                return densityWindows;
            }

            BrowseableBAMIndex browseableIndex = (BrowseableBAMIndex)index;

            // Calculate the first and last level 2 bin for this reference
            int firstBinNum = AbstractBAMFileIndex.getFirstBinInLevel(2);
            int lastBinNum = AbstractBAMFileIndex.getFirstBinInLevel(3) - 1;

            // For each level 2 bin (1Mb regions)
            for (int binNum = firstBinNum; binNum <= lastBinNum; binNum++) {
                Bin bin = new Bin(referenceIndex, binNum);

                // Get the genomic region this bin represents
                int binStart = browseableIndex.getFirstLocusInBin(bin);
                int binEnd = Math.min(browseableIndex.getLastLocusInBin(bin), sequenceLength);

                if (binStart >= sequenceLength) break;

                // Get overlapping region to find chunks
                BAMFileSpan span = index.getSpanOverlapping(referenceIndex, binStart, binEnd);
                if (span == null || span.getChunks() == null || span.getChunks().isEmpty()) continue;

                // Estimate reads based on chunk sizes
                long estimatedReads = 0;
                for (Chunk chunk : span.getChunks()) {
                    estimatedReads += (chunk.getChunkEnd() - chunk.getChunkStart());
                }

                if (estimatedReads > 0) {
                    densityWindows.add(new DensityWindow(binStart, binEnd, estimatedReads));
                }
            }

            return densityWindows;
        }
    }

    /**
     * Comparator for ReadEndsForMarkDuplicates that orders by read1 position then pair orientation then read2 position.
     */
    static class ReadEndsMDComparator implements Comparator<ReadEndsForMarkDuplicates> {

        final boolean useBarcodes;

        public ReadEndsMDComparator(final boolean useBarcodes) {
            this.useBarcodes = useBarcodes;
        }

        public int compare(final ReadEndsForMarkDuplicates lhs, final ReadEndsForMarkDuplicates rhs) {
            int compareDifference = lhs.libraryId - rhs.libraryId;
            if (useBarcodes) {
                final ReadEndsForMarkDuplicatesWithBarcodes lhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) lhs;
                final ReadEndsForMarkDuplicatesWithBarcodes rhsWithBarcodes = (ReadEndsForMarkDuplicatesWithBarcodes) rhs;
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.barcode, rhsWithBarcodes.barcode);
                }
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.readOneBarcode, rhsWithBarcodes.readOneBarcode);
                }
                if (compareDifference == 0) {
                    compareDifference = Integer.compare(lhsWithBarcodes.readTwoBarcode, rhsWithBarcodes.readTwoBarcode);
                }
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read1Coordinate - rhs.read1Coordinate;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.orientation - rhs.orientation;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
            }
            if (compareDifference == 0) {
                compareDifference = lhs.read2Coordinate - rhs.read2Coordinate;
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getTile() - rhs.getTile();
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getX() - rhs.getX();
            }

            if (compareDifference == 0) {
                compareDifference = lhs.getY() - rhs.getY();
            }

            // The following is arbitrary and is only included for completeness.
            // Other implementations may chose to forgo this tiebreak if they do not have
            // access to the index-in-file of the records (e.g. SPARK implmentations)

            if (compareDifference == 0) {
                compareDifference = Long.compare(lhs.read1IndexInFile, rhs.read1IndexInFile);
            }
            if (compareDifference == 0) {
                compareDifference = Long.compare(lhs.read2IndexInFile, rhs.read2IndexInFile);
            }

            return compareDifference;
        }
    }
}
