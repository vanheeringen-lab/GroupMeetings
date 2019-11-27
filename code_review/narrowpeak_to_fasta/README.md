# narrowpeak to fasta
### Introduction
I downloaded ATAC-seq data for many different species and developmental stages. I have aligned these
and called peaks on those. Hopefully the regulatory code between species is conserved and we can
infer for instance (conserved) transcription factor motifs, (conserved) chromatin remodelling, etc.  

To get a quick feeling for how conserved ATAC-seq peaks are between species, and to get a reference
benchmark we want to predict based purely on a sequence (of for instance 100bp) whether or not this
sequence belongs to a peak (and thus whether or not it is open). 

### Problem
I have peak files of different species (assemblies) and different stages, e.g.:
  - *mm10-early_2cell_peaks.narrowPeak*
  - *mm10-4cell_peaks.narrowPeak*
  - *mm10-8cell_peaks.narrowPeak*
  - *dm6-4hrs_peaks.narrowPeak*
  - *etc.*

I want to have a (sort-of) database where I store the sequence of each peak. I want to be able to 
quickly get all the peaks that belong to a sample, to an assembly, to a chromosome. 

Ideally I could dynamically get these peaks from a database (or data structure), and simply change 
a parameter (e.g. how large the sequence is) and that I do not have the regenerate everything from
scratch.

### Current solution
**peaks:**

I find each peak (black line) and take the sequence (green line) around the summit (red bar): 
<p align="center">
    <img src="https://raw.githubusercontent.com/vanheeringen-lab/GroupMeetings/master/code_review/narrowpeak_to_fasta/stagetofastq.jpg">
</p>

**non-peaks:**

I take the complement (red lines) of the peaks (black lines), and randomly sample the same amount of
sequences as there are peaks: 
<p align="center">
    <img src="https://raw.githubusercontent.com/vanheeringen-lab/GroupMeetings/master/code_review/narrowpeak_to_fasta/stagestofastq.jpg">
</p>

See the code [here](https://github.com/vanheeringen-lab/GroupMeetings/blob/master/code_review/narrowpeak_to_fasta/peak_to_fasta.py)
