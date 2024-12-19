package picard.sam.markduplicates.util;

import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;

import java.io.*;
import java.nio.BufferUnderflowException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.channels.WritableByteChannel;

/**
 * Optimized codec for ReadEnds that uses memory mapping for improved performance.
 * Falls back to standard I/O for non-file streams.
 */
public class ReadEndsForMarkDuplicatesCodec implements SortingCollection.Codec<ReadEndsForMarkDuplicates> {
    private static final int CHUNK_SIZE = 64 * 1024 * 1024; // 64MB chunks
    private static final int WRITE_BUFFER_SIZE = 8 * 1024 * 1024; // 8MB write buffer

    // Input handling
    private FileChannel inputChannel;
    private ByteBuffer mappedInputBuffer;
    private DataInputStream standardInput;
    private long currentInputPosition = 0;

    // Output handling
    private FileChannel outputChannel;
    private ByteBuffer writeBuffer;
    private DataOutputStream standardOutput;

    // Record size tracking for pre-allocation
    private static final int ESTIMATED_RECORD_SIZE = 50; // Approximate size in bytes

    @Override
    public SortingCollection.Codec<ReadEndsForMarkDuplicates> clone() {
        return new ReadEndsForMarkDuplicatesCodec();
    }

    @Override
    public void setOutputStream(final OutputStream os) {
        if (os instanceof FileOutputStream) {
            try {
                outputChannel = ((FileOutputStream) os).getChannel();
                writeBuffer = ByteBuffer.allocateDirect(WRITE_BUFFER_SIZE);
            } catch (Exception e) {
                // Fall back to standard output
                outputChannel = null;
                standardOutput = new DataOutputStream(new BufferedOutputStream(os, WRITE_BUFFER_SIZE));
            }
        } else {
            standardOutput = new DataOutputStream(new BufferedOutputStream(os, WRITE_BUFFER_SIZE));
        }
    }

    @Override
    public void setInputStream(final InputStream is) {
        if (is instanceof FileInputStream) {
            try {
                inputChannel = ((FileInputStream) is).getChannel();
                mapNextInputChunk();
            } catch (IOException e) {
                throw new PicardException("Failed to initialize memory-mapped input", e);
            }
        } else {
            standardInput = new DataInputStream(new BufferedInputStream(is, WRITE_BUFFER_SIZE));
        }
    }

    private void mapNextInputChunk() throws IOException {
        if (currentInputPosition >= inputChannel.size()) {
            mappedInputBuffer = null;
            return;
        }

        long remainingBytes = inputChannel.size() - currentInputPosition;
        long chunkSize = Math.min(CHUNK_SIZE, remainingBytes);

        mappedInputBuffer = inputChannel.map(
                FileChannel.MapMode.READ_ONLY,
                currentInputPosition,
                chunkSize
        );
        currentInputPosition += chunkSize;
    }

    private void ensureWriteBufferCapacity(int required) throws IOException {
        if (writeBuffer.remaining() < required) {
            writeBuffer.flip();
            outputChannel.write(writeBuffer);
            writeBuffer.clear();
        }
    }

    @Override
    public void encode(final ReadEndsForMarkDuplicates read) {
        try {
            if (outputChannel != null) {
                encodeToBuffer(read);
            } else {
                encodeToStream(read);
            }
        } catch (IOException ioe) {
            throw new PicardException("Exception writing ReadEnds to file.", ioe);
        }
    }

    private void encodeToBuffer(final ReadEndsForMarkDuplicates read) throws IOException {
        // Ensure we have enough space
        ensureWriteBufferCapacity(ESTIMATED_RECORD_SIZE);

        // Write all fields to buffer
        writeBuffer.putShort(read.score);
        writeBuffer.putShort(read.libraryId);
        writeBuffer.put(read.orientation);
        writeBuffer.putInt(read.read1ReferenceIndex);
        writeBuffer.putInt(read.read1Coordinate);
        writeBuffer.putLong(read.read1IndexInFile);
        writeBuffer.putInt(read.read2ReferenceIndex);

        if (read.orientation > ReadEnds.R) {
            writeBuffer.putInt(read.read2Coordinate);
            writeBuffer.putLong(read.read2IndexInFile);
        }

        writeBuffer.putShort(read.readGroup);
        writeBuffer.putShort(read.tile);
        writeBuffer.putShort((short)read.x);
        writeBuffer.putShort((short)read.y);
        writeBuffer.put(read.orientationForOpticalDuplicates);
        writeBuffer.putInt(read.duplicateSetSize);
        writeBuffer.put(read.firstOfFlag);
        writeBuffer.putInt(read.windowIndex);
    }

    private void encodeToStream(final ReadEndsForMarkDuplicates read) throws IOException {
        standardOutput.writeShort(read.score);
        standardOutput.writeShort(read.libraryId);
        standardOutput.writeByte(read.orientation);
        standardOutput.writeInt(read.read1ReferenceIndex);
        standardOutput.writeInt(read.read1Coordinate);
        standardOutput.writeLong(read.read1IndexInFile);
        standardOutput.writeInt(read.read2ReferenceIndex);

        if (read.orientation > ReadEnds.R) {
            standardOutput.writeInt(read.read2Coordinate);
            standardOutput.writeLong(read.read2IndexInFile);
        }

        standardOutput.writeShort(read.readGroup);
        standardOutput.writeShort(read.tile);
        standardOutput.writeShort((short)read.x);
        standardOutput.writeShort((short)read.y);
        standardOutput.writeByte(read.orientationForOpticalDuplicates);
        standardOutput.writeInt(read.duplicateSetSize);
        standardOutput.writeByte(read.firstOfFlag);
        standardOutput.writeInt(read.windowIndex);
    }

    @Override
    public ReadEndsForMarkDuplicates decode() {
        try {
            if (inputChannel != null) {
                return decodeFromBuffer();
            } else {
                return decodeFromStream();
            }
        } catch (EOFException eof) {
            return null;
        } catch (IOException ioe) {
            throw new PicardException("Exception reading ReadEnds from file.", ioe);
        }
    }

    private ReadEndsForMarkDuplicates decodeFromBuffer() throws IOException {
        // Check if we need to map next chunk
        if (mappedInputBuffer == null || !mappedInputBuffer.hasRemaining()) {
            mapNextInputChunk();
            if (mappedInputBuffer == null) {
                return null; // EOF
            }
        }

        try {
            final ReadEndsForMarkDuplicates read = new ReadEndsForMarkDuplicates();

            read.score = mappedInputBuffer.getShort();
            read.libraryId = mappedInputBuffer.getShort();
            read.orientation = mappedInputBuffer.get();
            read.read1ReferenceIndex = mappedInputBuffer.getInt();
            read.read1Coordinate = mappedInputBuffer.getInt();
            read.read1IndexInFile = mappedInputBuffer.getLong();
            read.read2ReferenceIndex = mappedInputBuffer.getInt();

            if (read.orientation > ReadEnds.R) {
                read.read2Coordinate = mappedInputBuffer.getInt();
                read.read2IndexInFile = mappedInputBuffer.getLong();
            }

            read.readGroup = mappedInputBuffer.getShort();
            read.tile = mappedInputBuffer.getShort();
            read.x = mappedInputBuffer.getShort();
            read.y = mappedInputBuffer.getShort();
            read.orientationForOpticalDuplicates = mappedInputBuffer.get();
            read.duplicateSetSize = mappedInputBuffer.getInt();
            read.firstOfFlag = mappedInputBuffer.get();
            read.windowIndex = mappedInputBuffer.getInt();

            return read;
        } catch (BufferUnderflowException e) {
            // If we hit the end of the buffer mid-record, map next chunk and retry
            mapNextInputChunk();
            if (mappedInputBuffer == null) {
                return null; // EOF
            }
            return decodeFromBuffer();
        }
    }

    private ReadEndsForMarkDuplicates decodeFromStream() throws IOException {
        final ReadEndsForMarkDuplicates read = new ReadEndsForMarkDuplicates();

        try {
            read.score = standardInput.readShort();
        } catch (EOFException eof) {
            return null;
        }

        read.libraryId = standardInput.readShort();
        read.orientation = standardInput.readByte();
        read.read1ReferenceIndex = standardInput.readInt();
        read.read1Coordinate = standardInput.readInt();
        read.read1IndexInFile = standardInput.readLong();
        read.read2ReferenceIndex = standardInput.readInt();

        if (read.orientation > ReadEnds.R) {
            read.read2Coordinate = standardInput.readInt();
            read.read2IndexInFile = standardInput.readLong();
        }

        read.readGroup = standardInput.readShort();
        read.tile = standardInput.readShort();
        read.x = standardInput.readShort();
        read.y = standardInput.readShort();
        read.orientationForOpticalDuplicates = standardInput.readByte();
        read.duplicateSetSize = standardInput.readInt();
        read.firstOfFlag = standardInput.readByte();
        read.windowIndex = standardInput.readInt();

        return read;
    }

    /**
     * Closes any open resources and ensures all data is written.
     */
    public void close() throws IOException {
        if (writeBuffer != null && writeBuffer.position() > 0) {
            writeBuffer.flip();
            outputChannel.write(writeBuffer);
        }

        if (standardOutput != null) {
            standardOutput.flush();
            standardOutput.close();
        }

        if (inputChannel != null) {
            inputChannel.close();
        }

        if (standardInput != null) {
            standardInput.close();
        }
    }
}