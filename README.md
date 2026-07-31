# INV_DUP long-read analysis pipeline

## 1. Project overview

This repository contains the pipeline I developed during my internship to analyse long-read alignments around known inversion regions and look for patterns compatible with inverted duplications (INV_DUP).

The analysis starts from inversion coordinates already defined in allcoords.txt and processes one BAM at a time. The main idea was to go step by step from the alignments found around each inversion to a more complete interpretation of the read structure.

I initially used the regional alignments to look for patterns that could distinguish INV_DUP reads. However, while inspecting the reads I found that some relevant supplementary alignments were missing from the regional fetch. Because of this, I added SA tag recovery before calculating the complete read geometry.

From there, the pipeline gradually became divided into several levels:

* first, candidate reads are selected from the alignments found around the inversion
* supplementary alignments are recovered when they provide information missing from the regional alignments
* all useful alignments are represented as ordered read segments
* the geometry of each read is calculated using strand, query gaps and overlaps, and reference coordinates
* each read receives a provisional classification
* POSSIBLE_INV_DUP reads are compared to see whether several reads support the same coordinate-level event
* finally, this event is compared with repeated regions already present in the reference, because some apparent INV_DUP patterns may already be explained by the reference sequence

I kept these steps separated because it made it easier to understand why a read was classified in a certain way and to change thresholds or rules without losing the original alignment information.

---

## 2. General pipeline

The main analysis follows these steps:

1. BAM and allcoords.txt are used to select candidate reads around the inversion.
2. SA tags are inspected and supplementary alignments are recovered when needed.
3. The alignments are converted into ordered read segments.
4. Read geometry and segment pairs are calculated.
5. Each read receives a provisional classification.
6. POSSIBLE_INV_DUP reads are compared to identify coordinate events supported by several reads.
7. A BAM-level summary is generated.
8. The main candidate event is compared with self-BLAST repeat pairs.

The self-BLAST analysis describes repeated sequence blocks already present around the inversion in the reference genome and is used in the final comparison with the read-supported event.

Each BAM is analysed independently. Several BAMs may correspond to the same inversion, so they use the same inversion reference region even if their processed self-BLAST files are organised by BAM.

---

## 3. Inputs and main parameters

The main inputs are:

* BAM/allcoords.txt, containing the inversion identifiers and coordinates
* REF/Std_ref2.fa, used as reference genome
* one BAM file to analyse

Shared paths and parameters are stored in scripts/config.py.

Some of the main parameters are:

* MARGIN = 50000
* MIN_IDENTITY = 90.0
* MIN_LENGTH = 1000
* MIN_GEOMETRY_MAPQ = 55
* MIN_SA_RECOVERY_MAPQ = 20
* MIN_GAP_COVERAGE_FRACTION = 0.80
* COORDINATE_MARGIN = 50
* MIN_READS_INV_DUP = 3

These values are still provisional and should be validated with more cases.

### Coordinate conventions

The BAM and BLAST parts do not originally use the same coordinate convention.

BAM-derived reference coordinates and corrected query intervals are treated as 0-based half-open intervals.

SA positions are 1-based, so SA_POS is converted to REF_START by subtracting 1 before the segment is stored.

Raw BLAST coordinates are 1-based inclusive. They are first converted from the local FASTA position back to genomic coordinates. Before comparing self-BLAST pairs with BAM events, the BLAST coordinates are converted to the 0-based half-open convention used in the BAM analysis.

BED files use the standard BED convention with a 0-based start.

---

## 4. Self-BLAST repeat detection

The self-BLAST part searches for repeated sequence blocks around each inversion.

Its purpose is not to classify INV_DUP reads directly. I use it later to check whether a candidate pattern detected from the reads could already be explained by repeated or inverted sequence present in the reference.

### Generate the inversion FASTA

run_self_blast.py reads the inversion coordinates and extracts the inversion plus MARGIN on both sides.

The extracted region goes from inversion start - MARGIN to inversion end + MARGIN.

If the start becomes smaller than 1, it is corrected to 1. The FASTA is stored in results/self_blast/FASTA/.

### Run and filter the self-BLAST

Each FASTA is aligned against itself with blastn.

The filtered hits must satisfy:

* identity >= MIN_IDENTITY
* alignment length >= MIN_LENGTH
* the hit cannot be a self-hit

The orientation is determined from the subject coordinates:

* direct if the subject coordinates go forward
* inverted if they go backwards

### Deduplicate reciprocal hits

The same repeated pair is normally reported twice by self-BLAST: once with block A as the query and block B as the subject, and once in the opposite direction.

I normalise and sort the two intervals so both hits produce the same pair key. Orientation is also included in the key, so a direct and an inverted pair are not mixed.

At the moment, this deduplication uses exact coordinates. Similar BLAST hits with slightly different boundaries may therefore remain as separate pairs.

### Repeat pairs and BED

Each deduplicated pair stores the two repeated blocks, orientation and a representative BLAST hit.

The local BLAST positions are converted back to genomic coordinates using the start of the extracted reference region.

generate_bed.py then produces two BED rows per pair, one for each block, which is useful for visual inspection.

The main processed outputs are stored under results/self_blast/OUTPUT_TSV/BAM_NAME/ and include:

* BAM_NAME_blast_hits.tsv
* BAM_NAME_blast_pairs.tsv
* BAM_NAME_blast_summary.tsv

---

## 5. Candidate reads and SA recovery

The first BAM analysis selects reads found around the inversion using their regional alignment patterns.

Patterns such as POSS_INV_DUP, SPLIT, SPLIT_INV and ONLY_SUPPL are only used to decide which reads continue to the detailed analysis. They are not the final INV_DUP classification.

I also kept ONLY_SUPPL reads because, while inspecting the BAMs, I found cases where useful reads were mainly represented by supplementary alignments and would otherwise be lost too early.

### Why I added SA recovery

At first I only worked with the alignments returned by the regional BAM fetch.

While checking individual reads, I saw that some SA tags described alignments that were not present in that regional set. This was especially important when two regional segments left a gap in the read and an additional supplementary alignment could explain what was inside that gap.

Because of this, the pipeline checks the SA tags before reconstructing the final segment geometry.

### SA processing

SA entries contain chromosome, position, strand, CIGAR and MAPQ.

Before trying to recover anything, the pipeline checks whether the SA alignment is already present among the regional alignments. Only genuinely new SA entries are considered for recovery.

New SA entries with MAPQ >= 20 are considered reliable enough to keep as supplementary evidence.

If the full alignment record can be found in the BAM, it is stored as RECOVERED_SA.

If the tag is reliable but the complete alignment record is not available in the BAM subset, the information is still stored as SA_TAG_ONLY.

This was important because the BAMs used during the internship were subsets of larger individual BAMs. A missing full record does not necessarily mean that the alignment never existed.

Each read also receives an SA recovery status such as NO_SA, SA_ALREADY_IN_REGION, NO_RELIABLE_SA, SA_RECOVERED or a status indicating that a reliable SA could not be recovered.

I kept the recovered SA information separate from the original regional pattern so that recovering a new alignment does not change what was initially observed in the inversion region.

---

## 6. Read segment representation

generate_read_segments.py converts all the useful alignment information into a common segment table.

Each row represents one alignment segment of one read. The SOURCE column indicates where the segment came from:

* REGION: already present in the regional BAM fetch
* RECOVERED_SA: recovered from the BAM using an SA tag
* SA_TAG_ONLY: reconstructed from the SA tag because the complete BAM record was not available

For SA_TAG_ONLY segments, some fields such as FLAG or raw pysam query coordinates are unknown. These empty values should not be interpreted as False.

### Corrected query coordinates

One of the main problems I had to solve was placing all alignments of the same read in the same query orientation.

Hard clipping complicates this because hard-clipped bases belong to the original read but are not stored in the BAM record sequence.

For forward alignments, the leading hard clipping is added to the query coordinates.

For reverse alignments, the interval is reflected and then shifted using the hard clipping that corresponds to the beginning of the original read.

The corrected coordinates are stored as CORR_QUERY_START and CORR_QUERY_END. These are the coordinates used later to order the segments and calculate gaps or overlaps.

All segments from the same read are sorted by corrected query coordinates and receive an ALIGNMENT_INDEX.

### Why low-MAPQ segments are still stored

I keep all candidate segments in read_segments.tsv, even when their MAPQ is too low for the main geometry.

The geometry filter is applied later. I preferred this because otherwise changing the MAPQ threshold would require regenerating the segment table, and useful supporting information could be lost too early.

---

## 7. Read geometry analysis

The geometry stage reconstructs how the read segments relate to each other.

I use a stricter MAPQ threshold for the alignments that define the main geometry: MIN_GEOMETRY_MAPQ = 55.

Reliable SA-derived alignments can still be used as supporting information with MIN_SA_RECOVERY_MAPQ = 20.

This separation was useful because an SA alignment may not be strong enough to define the main geometry but can still explain a gap between two high-confidence regional alignments.

### Segment pairs

The pipeline stores all possible pairs of geometry segments and also records whether each pair is consecutive.

I initially focused on consecutive segments, but later I found that an inverted segment could overlap in reference coordinates with a later alignment even when another segment was located between them in query order.

For this reason, reference-overlap compatibility can be checked between non-consecutive segments as well.

Query gaps and query overlaps, however, are calculated from consecutive segments because they describe continuity along the read.

The main pair information includes:

* whether the segments are on the same chromosome
* whether they have opposite strands
* query gap or overlap
* reference overlap
* reference overlap fraction
* reference length ratio
* whether they are consecutive

### Query gaps

A query gap is not automatically interpreted as missing sequence.

While inspecting reads, I found that an apparent gap could be explained by supplementary alignments that were not part of the main regional geometry.

Because of this, the pipeline checks whether reliable SA-derived segments cover the original gap between two REGION alignments.

### Multiple SA alignments: HsInv1146

At first I looked for one supplementary alignment that explained the gap.

This did not work well for HsInv1146. In this case, some gaps were explained by several supplementary alignments together.

I changed the logic so that the query coverage of all reliable supplementary segments inside the REGION gap is merged. The gap is considered supported when their combined coverage reaches at least 80% of the original gap and at least one of the supporting segments is inverted relative to the REGION flanks.

This allows patterns such as:

* REGION + | SA - | REGION +
* REGION + | SA + | SA - | REGION +
* REGION + | SA - | SA - | REGION +

The supporting sequence may come from the same chromosome or another chromosome, so I keep intrachromosomal and interchromosomal gap cases as separate classification reasons.

---

## 8. Read classification

classify_inv.py assigns each read to one of four provisional groups:

* POSSIBLE_INV_DUP
* NOT_INV_DUP
* AMBIGUOUS
* NOT_ENOUGH_INFORMATION

I call this classification provisional because the read still has to be compared with other reads at coordinate level.

Some of the main rules are:

### Not enough information

A read with fewer than two reliable geometry segments is classified as NOT_ENOUGH_INFORMATION and INSUFFICIENT_GEOMETRY.

One segment is simply not enough to decide.

### Large query overlap

A large query overlap is currently classified as NOT_INV_DUP and LARGE_QUERY_OVERLAP.

I interpreted these cases as more compatible with redundant or ambiguous alignments than with independent structural segments. I also observed this type of pattern in some of the normal inversion controls.

### Same-strand reads

Same-strand geometry alone does not support an inverted structure.

However, if an inverted SA-derived alignment explains the gap between the REGION segments, the read can still become POSSIBLE_INV_DUP and REGION_GAP_WITH_INVERTED_SA_SUPPORT.

A large same-strand gap without this support is kept as AMBIGUOUS rather than automatically rejected.

### Mirror-like reads

Strongly symmetric opposite-strand pairs with very high reference overlap and similar lengths are classified as AMBIGUOUS and MIRROR_LIKE_GEOMETRY.

I kept them as ambiguous because these patterns may represent redundant or uncertain mappings.

### Local overlap

A read with opposite-strand segments on the same chromosome and reference overlap can be classified as POSSIBLE_INV_DUP and COMPATIBLE_INTRACHROMOSOMAL_OVERLAP.

For exactly two segments, the equivalent classification is POSSIBLE_INV_DUP and COMPATIBLE_TWO_SEGMENT_PATTERN.

If the two inverted segments do not overlap in reference coordinates, the current two-segment model does not consider them compatible.

### SA-supported gaps

For reads with more complex geometry, an inverted supplementary structure explaining the original query gap can produce POSSIBLE_INV_DUP and DISTANT_INTRACHROMOSOMAL_GAP_WITH_INVERTED_SA or POSSIBLE_INV_DUP and INTERCHROMOSOMAL_GAP_WITH_INVERTED_SA.

Mixed-strand cases without enough supporting evidence are kept as AMBIGUOUS.

---

## 9. Coordinate events across reads

Read-level classification was not enough.

While testing the rules, I saw that even normal inversion controls could contain individual reads classified as POSSIBLE_INV_DUP. Because of this, I added another condition: several candidate reads should point to approximately the same genomic event.

classify_coords.py takes only POSSIBLE_INV_DUP reads and builds coordinate signatures.

### Local overlap events

For local overlap candidates, the compatible pair with the largest reference overlap is selected and its overlap coordinates are stored.

If several reads have similar coordinates, they are grouped into the same event.

I found that the complete overlap interval was not always equally stable between reads. In HsInv0231, for example, one boundary was much more recurrent than the complete overlap length.

Because of this, the pipeline compares the start and end recurrence and uses the more stable boundary as the main event anchor.

COORDINATE_MARGIN is currently 50 bp.

### GAP events

For SA-supported gap candidates, the signature includes the REGION boundaries and the SA support.

Two GAP reads are grouped only if the regional coordinate is compatible and their SA support is also similar.

Multiple SA alignments are kept because complex reads such as HsInv1146 cannot always be represented by one supplementary alignment.

An event is considered recurrent when at least two reads support it.

### Representative event coordinates
The reads grouped into the same event do not always have exactly the same coordinates. Small differences can appear between alignments, so I do not use the coordinates of a single read as the final event position.
For each type of event, the pipeline first checks whether the start or the end coordinate is more recurrent between the candidate reads. Two coordinates are considered compatible when they are within COORDINATE_MARGIN. The most recurrent boundary is then used as the anchor for grouping the reads.
As reads are added to an event, the event anchor is updated using the median of the anchor coordinates from all the reads already assigned to that event. I used the median because it gives a simple representative coordinate while being less affected by a read with a slightly more distant breakpoint.
Once the events have been grouped, the final coordinate table stores the median start and median end across all reads supporting each event. Therefore, these values represent the central coordinates of the event and do not necessarily correspond exactly to the interval of one individual read.
For GAP events, the same idea is also applied to the supplementary alignments. MEDIAN_SA_START and MEDIAN_SA_END summarize the coordinates of the SA segments supporting the event.
For now, I only keep these median coordinates in the event summary. This gives a simple representative position, but it does not show how dispersed the individual read coordinates are. A future improvement could also report the coordinate distribution, for example the 25th and 75th percentiles, to show the variability between supporting reads.

---

## 10. Final BAM summary and self-BLAST comparison

After the coordinate analysis, the pipeline creates one summary row for the BAM.

The main event is the coordinate event supported by the largest number of reads.

The summary includes, among other values:

* number of analysed reads
* number of POSSIBLE_INV_DUP reads
* number of coordinate events
* number of recurrent events
* main event type
* main event coordinates
* number of reads supporting the main event
* IS_POSS_INV_DUP

### MAIN_FRACTION_POSSIBLE_READS_IN_MAIN_EVENT

I added MAIN_FRACTION_POSSIBLE_READS_IN_MAIN_EVENT to see whether the POSSIBLE_INV_DUP reads were really describing the same event.

It is calculated as the number of reads supporting the main event divided by the total number of POSSIBLE_INV_DUP reads.

For example, 15 reads in the main event out of 20 POSSIBLE_INV_DUP reads gives 0.75.

This is only a simple measure of consistency between reads, not a statistical score.

### BAM-level candidate call

The current final call is:

* YES when there is a recurrent event and the BAM contains at least MIN_READS_INV_DUP reads classified as POSSIBLE_INV_DUP
* NOT_ENOUGH_READS when there is a recurrent event but the BAM does not contain enough POSSIBLE_INV_DUP reads
* NO when there is no recurrent coordinate event

MIN_READS_INV_DUP is currently 3.

### Final comparison with the reference

The self-BLAST is used only after the read-supported event has already been defined.

The goal is to detect cases where part of the apparent INV_DUP pattern could already be explained by an inverted repeat in the reference.

The BLAST pair coordinates are converted to the same coordinate convention used for the BAM event before the comparison.

At the moment, I keep this step conservative.

A nearby reference repeat does not automatically mean that the candidate is false. Therefore, the current result is either KEEP or REVIEW_REFERENCE_REPEAT rather than automatically changing the candidate to NOT_INV_DUP.

---

## 11. Outputs

For each BAM, the read-analysis results are stored in results/invdup_filtering/BAM_NAME/.

The most important files are:

* BAM_NAME_read_segments.tsv: one row per alignment segment
* BAM_NAME_segments_geometry.tsv: one geometry summary per read
* BAM_NAME_pairs.tsv: pairwise geometry between reliable segments
* BAM_NAME_classification.tsv: provisional classification for every read
* BAM_NAME_coordinate_events.tsv: recurrent coordinate events from POSSIBLE_INV_DUP reads
* BAM_NAME_summary.tsv: final BAM-level summary

The summary printed in the terminal only shows the most important values. The TSV files should be checked when the supporting READ_IDs or intermediate information are needed.

The self-BLAST outputs are stored under results/self_blast/OUTPUT_TSV/BAM_NAME/.

---

## 12. Running the pipeline

The final pipeline processes one BAM at a time.

From the project root:

python3 -m scripts.run_pipeline HsInv0231_2.bam

The BAM is searched inside BAM_DIR and the corresponding output folder is created automatically.

The self-BLAST analysis generates the repeat pairs used for the final comparison with the read-supported event.

---

## 13. Environment

The main Python dependencies are pandas and pysam.

A minimal environment can be created with:

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

---

## 14. Things that could still be improved

The pipeline is now complete enough to run from a BAM to a final summary, but some decisions are still provisional.

* The MAPQ, query gap, query overlap, coordinate margin and minimum read thresholds should be validated with more cases.
* The different MAPQ thresholds for main geometry and SA supporting evidence should be confirmed.
* MIN_READS_INV_DUP should probably be reviewed to decide whether it should use all POSSIBLE_INV_DUP reads or only the reads supporting the main event.
* Self-BLAST pairs are currently deduplicated using exact coordinates, so very similar hits can remain as different pairs.
* The final reference-repeat comparison currently marks suspicious cases for review instead of discarding them automatically.
* GAP events with several SA alignments can be structurally complex, and median SA coordinates may not describe every case well.

*Initial development by Núria Regué Sorribes during an internship at the Comparative and Functional Genomics group (IMIM).*
