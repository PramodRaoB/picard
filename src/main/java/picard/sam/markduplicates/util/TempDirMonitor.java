package picard.sam.markduplicates.util;

import java.io.File;
import java.nio.file.*;
import java.io.IOException;
import htsjdk.samtools.util.Log;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

public class TempDirMonitor {
    private static final Log log = Log.getInstance(TempDirMonitor.class);
    private final List<File> tempDirs;

    public TempDirMonitor(List<File> tempDirs) throws IOException {
        // Create a dedicated temp directory
        this.tempDirs = tempDirs;
        log.debug("Monitoring temp directories");
    }

    public void logCurrentState() {
        try {
            log.info("\n=== Current Temporary Directory State ===");
            Map<Path, Long> fileSizes = new HashMap<>();
            AtomicInteger cnt = new AtomicInteger();

            // Collect all file information
            for (File tempDir : tempDirs) {
                Files.walk(tempDir.toPath())
                        .filter(Files::isRegularFile)
                        .forEach(p -> {
                            try {
                                cnt.getAndIncrement();
                                fileSizes.put(p, Files.size(p));
                            } catch (IOException e) {
                                log.warn("Could not get size for: " + p);
                            }
                        });
            }

            // Log total size
            long totalSize = fileSizes.values().stream().mapToLong(Long::valueOf).sum();
            log.info(String.format("Total directory size: %s (%d bytes) (%d files)",
                    humanReadableByteCount(totalSize), totalSize, cnt.get()));

            log.info("\n=====================================");

        } catch (IOException e) {
            log.warn("Could not measure directory state", e);
        }
    }

    private String humanReadableByteCount(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String pre = "KMGTPE".charAt(exp-1) + "";
        return String.format("%.1f %sB", bytes / Math.pow(1024, exp), pre);
    }
}