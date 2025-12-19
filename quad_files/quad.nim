# Andreas Wilm <andreas.wilm@gmail.com>
# Copyright 2023
# License: MIT


# Reads a BAM file maps base qualities as defined in a base quality
# mapping file (CSV)


import std/strformat
import std/tables
import std/strutils

import csvtools
import hts


proc parse_qmap_file(qmap_file: string): Table[uint8, uint8] =
  # parse the quality mapping csv file and turn into a table.
  # file is expected to have a header.
  # values should be lsited as: qold,qnew

  let max_val = 93
  var qmap = initTable[uint8, uint8]()
  type Q2Q = object
    key: int
    value: int

  for i in countup(0, int(high(uint8))): # nim requires the int cast, but why?
    qmap[uint8(i)] = uint8(i)


  for row in csv[Q2Q](qmap_file, skip_header = true):
    for v in [row.key, row.value]:
      doAssert(v >= 0 and v <= max_val, fmt"Invalid quality {v} found in {qmap_file}")
    qmap[uint8(row.key)] = uint8(row.value)
  return qmap


proc main(qmap_file: string, in_bam: string, out_bam: string, out_fmt = "",
    threads = 1): int =
  var
    ibam: Bam
    obam: Bam

  let qmap = parse_qmap_file(qmap_file)

  if not open(ibam, in_bam, index = false):
    # htsnim prints an error
    return 1

  open(obam, out_bam, threads = threads, mode = "w" & out_fmt)
  obam.write_header(ibam.hdr)

  for aln in ibam:
    var bqs = bam_get_qual(aln.b)
    for i in 0..<aln.b.core.l_qseq:
      bqs[i] = qmap[bqs[i]]

    obam.write(aln)

  obam.close()
  return 0


when isMainModule:
  import cligen
  dispatch(main)
