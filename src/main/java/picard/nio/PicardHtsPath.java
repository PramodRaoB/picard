/*
 * The MIT License
 *
 * Copyright (c) 2022 The Broad Institute
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

package picard.nio;

import htsjdk.io.HtsPath;
import htsjdk.io.IOPath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.utils.ValidationUtils;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Collectors;

/**
 * A Subclass of {@link HtsPath} with conversion to {@link Path} making use of {@link IOUtil}
 */
public class PicardHtsPath extends HtsPath {
    /**
     * Create a PicardHtsPath from a raw input path string.
     * <p>
     * If the raw input string already contains a scheme (including a "file" scheme), assume its already
     * properly escape/encoded. If no scheme component is present, assume it references a raw path on the
     * local file system, so try to get a Path first, and then retrieve the URI from the resulting Path.
     * This ensures that input strings that are local file references without a scheme component and contain
     * embedded characters are valid in file names, but which would otherwise be interpreted as excluded
     * URI characters (such as the URI fragment delimiter "#") are properly escape/encoded.
     *
     * @param rawInputString a string specifying an input path. May not be null.
     */
    public PicardHtsPath(final String rawInputString) {
        super(rawInputString);
    }

    /**
     * Create a PicardHtsPath from an existing {@link HtsPath} or subclass.
     *
     * @param htsPath an existing PathSpecifier. May not be null.
     */
    public PicardHtsPath(final HtsPath htsPath) {
        super(htsPath);
    }

    /**
     * Create a PicardHtsPath from a {@link File} reference. Uses the URI string of {@code file}.
     * @param file the file reference to create this object from
     */
    public PicardHtsPath(final File file){
        this(file.toURI().toString());
    }

    /**
     * Create a {@link List<PicardHtsPath>} from path representations.
     * @param paths URIs or local paths. May not be null but may be empty.
     * @return the converted {@link List}
     */
    public static List<PicardHtsPath> fromPaths(Collection<String> paths) {
        Objects.requireNonNull(paths);
        return paths.stream().map(PicardHtsPath::new).collect(Collectors.toList());
    }

    /**
     * Resolve the URI of this object to a {@link Path} object.
     *
     * @return the resulting {@link Path}
     * @throws RuntimeException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        try {
            return IOUtil.getPath(super.getURIString());
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Construct a {@link PicardHtsPath} from a {@link Path}
     * @param path may NOT be null
     * @return a new object representing path
     */
    public static PicardHtsPath fromPath(final Path path){
        Objects.requireNonNull(path);
        return new PicardHtsPath(path.toUri().toString());
    }

    /**
     * Test if {@code ioPath} is something other than a regular file, directory, or symbolic link.
     *
     * @return {@code true} if it's a device, named pipe, htsget API URL, etc.
     * @throws RuntimeException if an I/O error occurs when creating the file system
     */
    public static boolean isOther(final IOPath ioPath) {
        if(ioPath.isPath()) {
            try {
                return Files.readAttributes(ioPath.toPath(), BasicFileAttributes.class).isOther();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        } else {
            return true;
        }
    }

    /**
     * Create a {@link List<Path>} from {@link PicardHtsPath}s
     * @param picardHtsPaths may NOT be null
     * @return Path representations of the input picardHtsPaths
     */
    public static List<Path> toPaths(final Collection<PicardHtsPath> picardHtsPaths){
        Objects.requireNonNull(picardHtsPaths);
        return picardHtsPaths.stream().map(PicardHtsPath::toPath).collect(Collectors.toList());
    }

    /**
     * Takes an IOPath and returns a new PicardHtsPath object that keeps the same basename as the original but has
     * a new extension. If append is set to false, only the last component of an extension will be replaced.
     * e.g. "my.fasta.gz" -> "my.fasta.tmp"
     *
     * If the input IOPath was created from a rawInputString that specifies a relative local path, the new path will
     * have a rawInputString that specifies an absolute path. (This perhaps unwanted conversion occurs when we call PicardHtsPath.toPath().)
     *
     * Examples:
     *     - (test_na12878.bam, .bai) -> test_na12878.bai (append = false)
     *     - (test_na12878.bam, .md5) -> test_na12878.bam.md5 (append = true)
     *
     * @param path The original path
     * @param append If set to true, append the new extension to the original basename. If false, replace the original extension
     *               with the new extension. If append = false and the original name has no extension, an exception will be thrown.
     * @param newExtension A new file extension. Must include the leading dot e.g. ".txt", ".bam"
     * @return A new PicardHtsPath object with the new extension
     */
    public static PicardHtsPath replaceExtension(final IOPath path, final String newExtension, final boolean append){
        ValidationUtils.validateArg(newExtension.startsWith("."), "newExtension must start with a dot '.'");

        final String oldFileName = path.toPath().getFileName().toString();

        String newFileName;
        if (append){
            newFileName = oldFileName + newExtension;
        } else {
            final Optional<String> oldExtension = path.getExtension();
            if (oldExtension.isEmpty()){
                throw new PicardException("The original path must have an extension when append = false: " + path.getURIString());
            }
            newFileName = oldFileName.replaceAll(oldExtension.get() + "$", newExtension);
        }

        return PicardHtsPath.fromPath(path.toPath().resolveSibling(newFileName));
    }

    /**
     * Wrapper for Path.resolve()
     */
    public static PicardHtsPath resolve(final PicardHtsPath absPath, final String relativePath){
        return PicardHtsPath.fromPath(absPath.toPath().resolve(relativePath));
    }
}
