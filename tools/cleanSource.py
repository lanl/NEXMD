#!/usr/bin/python
#
# Call this with a list of files to consider.

import optparse, os, sys, re, tempfile

parser = optparse.OptionParser(description = """This script cleans a source
file. It does (1) remove trailing white space characters, and (2) replace a
TAB character with space characters in cases where the TAB character follows a
leading space character.""")

parser.add_option("--pretend",
    help = "Just check, but not change any files",
    action = "store_true",
    default = False,
    dest = "pretend")

parser.add_option("--tab",
    help = "Replace a TAB character with N spaces (default N = 8)",
    default = 8,
    type = "int",
    dest = "tab",
    metavar = "N")

options, arguments = parser.parse_args()

if len(arguments) == 0:
  print("what files should I fix?")
  sys.exit(1)

for file in arguments:
  if not os.path.isfile(file):
    continue

  try:
    fd = open(file)
    lines = fd.readlines()
    fd.close()

  except:
    print("error opening file " + file)
    sys.exit(1)

  lineNumber = 0
  fileNeedsFixing = False
  for lineNumber in range(len(lines)):
    if re.compile("[ \t]+$").search(lines[lineNumber]):
      fileNeedsFixing = True
      print("trailing whitespace in file " + file + " on line " + str(lineNumber+1) + ": " + lines[lineNumber].rstrip())
      if not options.pretend:
        lines[lineNumber] = lines[lineNumber].rstrip()

    if re.compile("^ +\t+").search(lines[lineNumber]):
      fileNeedsFixing = True
      print("indent SP followed by TAB in file " + file + " on line " + str(lineNumber+1) + ": " + lines[lineNumber].rstrip())

      lineCopy = lines[lineNumber]
      tabInSpace = ""
      for i in range(options.tab):
        tabInSpace += " "

      while True:
        result = re.compile("^( +)\t").search(lineCopy)
        if result:
          lineCopy = re.sub("^ +\t", result.group(1) + tabInSpace, lineCopy)
        else:
          break

      print("fixed SP followed by TAB in file  " + file + " on line " + str(lineNumber+1) + ": " + lineCopy.rstrip())

      if not options.pretend:
        line[lineNumber] = lineCopy

  if not options.pretend and fileNeedsFixing:
    print("writing out fixed file " + file)
    fd = open(file, "w")
    for line in lines:
      print >> fd, line.rstrip()
    fd.close()
