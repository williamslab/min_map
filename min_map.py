#!/usr/bin/env python3

"""
Given an input genetic map and a maximum error rate, greedily drops positions
to produce a map with fewer entries. Following processing, the only remaining
map positions are those that, if dropped, would produce an error greater than
the threshold (assuming linear interpolation for finding genetic positions).

Author: Amy Williams
Date: 21 Oct 2020

This program is distributed under the terms of the GNU General Public License
"""

import argparse


class MapEntry:
  def __init__(self, physPos, genetPos, lastEntry, genetPos2 = None):
    self.physPos = physPos
    self.genetPos = genetPos
    self.genetPos2 = genetPos2
    self.left = lastEntry
    self.right = None
    # Nodes to the right of <self> that have been deleted
    self.rightDel = set()
    self.error = 10000

def main():
  """ Handle command line options """
  parser = argparse.ArgumentParser(
      description='Minimum viable genetic map: deletes sites allowing up to -error tolerance')
  parser.add_argument('-error', action='store', dest='error', type=float, default=0.01, help='maximum allowable error in interpolating any map position (default 0.01 map units)')
  parser.add_argument('-chr', action='store', dest='chr', type=str, required=True, help='chromosome to analyze')
  parser.add_argument('-mapfile', action='store', dest='mapfile', type=str, required=True, help='map filename')
  parser.add_argument('-out', action='store', dest='prefix', type=str, default='min_viable_map', help='output prefix (default min_viable_map)')
  parser.add_argument('-physcol', action='store', dest='physcol', type=int, default=1, help='0-based column number of physical positons (default 1)')
  parser.add_argument('-genetcol', action='store', dest='genetcol', type=int, default=3, help='0-based column number of genetic position (default 3)')
  parser.add_argument('-noheader', action='store_true', dest='keepheader', default=False, help='input map is all entries, so keep the first line (default False)')
  # for the other sex's genetic position:
  parser.add_argument('-sexspecific', action='store_true', dest='sexspec', default=False, help='input map is sex-specific (default False)')
  parser.add_argument('-genet2col', action='store', dest='genet2col', type=int, help='used with -sexpsecific: 0-based column number of second sex\'s genetic position (no default)')

  args = parser.parse_args()

  if args.sexspec:
    if args.genet2col == None:
      print("ERROR: when using -sexspecific, must set -genetcol and -genet2col to the")
      print("       0-based column of the male/female positions")
      quit()
    elif args.genetcol == args.genet2col:
      print("ERROR: need to set -genetcol and/or -genet2col to _different_ 0-based columns")
      print("       of the male/female positions")
      quit()
    if args.genetcol > args.genet2col:
      print("ERROR: set -genetcol to an index below -genet2col (this ensures the output")
      print("       uses the same column order as the input)")
      quit()

  """ Read in genetic map """
  the_map = None
  lastEntry = None
  entryCount = 0

  with open(args.mapfile, "r") as f:
    if not args.keepheader:
      f.readline() # throw away header

    for line in f:
      if line[0] == '#':
        # comment
        continue

      fields = line.strip().split()
      physPos = int(fields[ args.physcol ])
      genetPos = float(fields[ args.genetcol ])
      if not args.sexspec: # sex averaged map:
        entry = MapEntry(physPos, genetPos, lastEntry)
      else: # sex-specific map:
        genetPos2 = float(fields[ args.genet2col ])
        entry = MapEntry(physPos, genetPos, lastEntry, genetPos2)

      if lastEntry != None:
        lastEntry.right = entry
      lastEntry = entry

      if the_map == None:
        the_map = entry

      entryCount += 1

  print("Pre-filter:  chrom {}: {} entries".format(args.chr, entryCount))


  """ Greedily drop positions until no more positions can be dropped """
  cur_entry = the_map
  while cur_entry != None:
    update_error(cur_entry)
    cur_entry = cur_entry.right

  # decide which entries to delete
  cur_entry = the_map
  first_skipped = None
  while cur_entry != None:
    if cur_entry.error >= args.error:
      # error above threshold: keep in map
      cur_entry = cur_entry.right
      continue

    # should only delete an entry if it has lower error than both its
    # neighbors
    shouldDel = (cur_entry.left == None or \
                  cur_entry.error <= cur_entry.left.error) and \
                (cur_entry.right == None or \
                  cur_entry.error <= cur_entry.right.error)

    if not shouldDel:
      # not going to delete <cur_entry> now -- one of its neighbors has lower
      # error
      if first_skipped == None:
        # might be deletable later (will delete others that have lower error
        # first), so revisit:
        first_skipped = cur_entry
      cur_entry = cur_entry.right
      continue

    # candidate for deletion: ensure that deleting <cur_entry> won't induce
    # an error greater than args.error in entries that were already deleted
    for prevDel in (cur_entry.left.rightDel | cur_entry.rightDel):
      if error(cur_entry.left, prevDel, cur_entry.right) >= args.error:
        shouldDel = False
        break

    if shouldDel:
      # <cur_entry> passes all checks: delete this entry
      cur_entry.left.right = cur_entry.right
      cur_entry.right.left = cur_entry.left
      update_error(cur_entry.left)
      update_error(cur_entry.right)

      cur_entry.left.rightDel |= cur_entry.rightDel
      cur_entry.left.rightDel.add(cur_entry)

      entryCount -= 1

      # move to next entry:
      if first_skipped != None:
        cur_entry = first_skipped
        first_skipped = None
      elif cur_entry.left.error < args.error:
        # corner case: left error can drop given the removal of <cur_entry>
        cur_entry = cur_entry.left
      else:
        cur_entry = cur_entry.right
    else:
      cur_entry = cur_entry.right

  print("Post-filter: chrom {}: {} entries".format(args.chr, entryCount))

  assert(the_map != None)

#  cur_entry = the_map
#  while cur_entry != None:
#    if cur_entry.error < 10000:
#      print("error: {}".format(cur_entry.error))
#    cur_entry = cur_entry.right

#  cur_entry = the_map
#  while cur_entry != None:
#    for prevDel in cur_entry.rightDel:
#      print("del-error: {}".format(error(cur_entry, prevDel, cur_entry.right)))
#    if cur_entry.right == None:
#      assert(cur_entry.rightDel == set())
#    cur_entry = cur_entry.right

  """ Print minimal map """
  outfile = "{}{}.txt".format(args.prefix, args.chr)
  with open(outfile, "w") as out:
    if not args.sexspec:
      print("# Chrom\tPosition\tRate\tMap(cM)", file=out)
    else:
      print("#chr\tpos\tmale_cM\tfemale_cM", file=out)
    cur_entry = the_map
    while cur_entry != None:
      if not args.sexspec:
        print("{}\t{}\tNA\t{}".format(args.chr, cur_entry.physPos,
                                      cur_entry.genetPos), file=out)
      else:
        print("{}\t{}\t{}\t{}".format(args.chr, cur_entry.physPos,
                                      cur_entry.genetPos, cur_entry.genetPos2),
                                      file=out)
      cur_entry = cur_entry.right



""" Updates the error if ```entry``` were interpolated from its current left
    and right positions
Arguments:
  entry: an entry in the genetic map
"""
def update_error(entry):
  if entry.left == None or entry.right == None:
    return

  entry.error = error(entry.left, entry, entry.right)

""" Calculates the error if ```target``` were interpolated from the ```left```
    and ```right``` entries
Arguments:
  left: an entry to the left of ```target`` in the genetic map
  target: an entry in the genetic map
  right: an entry to the right of ```target``` in the genetic map
"""
def error(left, target, right):
  # interpolate and calculate error for target position
  interpFrac = (target.physPos - left.physPos) / \
               (right.physPos  - left.physPos)
  interpGenet = left.genetPos + \
                interpFrac * (right.genetPos - left.genetPos)

  if left.genetPos2 == None:
    assert(right.genetPos2 == None) # sex avg map
    # error:
    return abs(target.genetPos - interpGenet)
  
  assert(right.genetPos2 != None)

  # sex specific map
  interpGenet2 = left.genetPos2 + \
                interpFrac * (right.genetPos2 - left.genetPos2)

  # return max error
  return max(abs(target.genetPos  - interpGenet),
             abs(target.genetPos2 - interpGenet2))


if __name__ == '__main__':
  main()
