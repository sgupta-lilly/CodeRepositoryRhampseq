import pysam

def is_low_error_probability(qualities, threshold=20):
    return all(q >= threshold for q in qualities)

def merge_overlapping_reads(bam_file, output_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    output_bam = pysam.AlignmentFile(output_file, "wb", template=bam)

    read_pairs = {}
    
    for read in bam.fetch(until_eof=True):
        qname = read.query_name
        if qname not in read_pairs:
            read_pairs[qname] = [None, None]
        if read.is_read1:
            read_pairs[qname][0] = read
        elif read.is_read2:
            read_pairs[qname][1] = read

    for read1, read2 in read_pairs.values():
        if read1 and read2:
            if read1.reference_start <= read2.reference_end and read2.reference_start <= read1.reference_end:
                overlap_start = max(read1.reference_start, read2.reference_start)
                overlap_end = min(read1.reference_end, read2.reference_end)
                overlap_length = overlap_end - overlap_start + 1

                if overlap_length > 5:
                    merged_seq = read1.query_sequence[:overlap_start - read1.reference_start]
                    merged_qual = read1.query_qualities[:overlap_start - read1.reference_start]

                    for i in range(overlap_start - read1.reference_start, overlap_end - read1.reference_start + 1):
                        base1 = read1.query_sequence[i]
                        base2 = read2.query_sequence[i - (overlap_start - read2.reference_start)]
                        qual1 = read1.query_qualities[i]
                        qual2 = read2.query_qualities[i - (overlap_start - read2.reference_start)]
                        
                        if qual1 >= qual2:
                            merged_seq += base1
                            merged_qual += qual1
                        else:
                            merged_seq += base2
                            merged_qual += qual2

                    merged_seq += read2.query_sequence[overlap_end - read2.reference_start + 1:]
                    merged_qual += read2.query_qualities[overlap_end - read2.reference_start + 1:]

                    read1.query_sequence = merged_seq
                    read1.query_qualities = merged_qual
                    output_bam.write(read1)
                else:
                    output_bam.write(read1)
                    output_bam.write(read2)
            else:
                output_bam.write(read1)
                output_bam.write(read2)
        elif read1:
            output_bam.write(read1)
        elif read2:
            output_bam.write(read2)

    bam.close()
    output_bam.close()

# Example usage
merge_overlapping_reads("input.bam", "output.bam")