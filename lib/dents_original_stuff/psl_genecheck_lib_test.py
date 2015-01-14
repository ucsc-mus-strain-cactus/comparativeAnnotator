"""
Test the sequence_lib, psl_lib classes and functions
Original author: Dent Earl
Modified by: Ian Fiddes
Ian Fiddes's own tests are in sequence_lib_tests.py
"""
from glob import glob
import os
import shutil
import string
import subprocess
import sys
import unittest
import lib.sequence_lib as seq_lib
import lib.psl_lib as psl_lib

def makeTempDirParent():
  """ make the parent temp dir directory
  """
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir'))


def removeTempDirParent():
  """ remove the parent temp dir directory
  """
  if os.path.exists(os.path.join(os.curdir, '.tempTestDir')):
    shutil.rmtree(os.path.join(os.curdir, '.tempTestDir'))


def makeTempDir(name=None):
  """ make the directory where all temporary test files will be stored.
  """
  makeTempDirParent()
  charSet = string.ascii_lowercase + '123456789'
  if name is None:
    while True:
      name = '%s_%s' % (''.join(random.choice(charSet) for x in xrange(4)),
                        ''.join(random.choice(charSet) for x in xrange(4)))
      if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
        break
  if not os.path.exists(os.path.join(os.curdir, '.tempTestDir', name)):
    os.mkdir(os.path.join(os.curdir, '.tempTestDir', name))
  return os.path.join(os.curdir, '.tempTestDir', name)


def removeDir(dirpath):
  """ destroy a directory
  """
  if os.path.exists(dirpath):
    shutil.rmtree(dirpath)
  if glob(os.path.join(os.path.dirname(dirpath), '*')) == []:
    # if this is the last tempDir to be destroyed, destroy the parent.
    removeTempDirParent()


def which(program):
  """which() acts like the unix utility which, but is portable between os.
  If the program does not exist in the PATH then 'None' is returned.
  """
  def is_exe(fpath):
    return os.path.exists(fpath) and os.access(fpath, os.X_OK)

  fpath, fname = os.path.split(program)
  if fpath != '':
    if is_exe(program):
      return program
  else:
    for path in os.environ['PATH'].split(os.pathsep):
      exe_file = os.path.join(path, program)
      if is_exe(exe_file):
        return exe_file
  return None


def runCommands(cmds, localTempDir, inPipes=None, outPipes=None, errPipes=None):
  """ Run commands from CMDS list.
  """
  if inPipes is None:
    inPipes = [None] * len(cmds)
  if outPipes is None:
    outPipes = [None] * len(cmds)
  if errPipes is None:
    errPipes = [None] * len(cmds)
  for i, c in enumerate(cmds, 0):
    if inPipes[i] is None:
      sin = None
    else:
      sin = subprocess.PIPE
    if outPipes[i] is None:
      sout = None
    else:
      sout = subprocess.PIPE
    if errPipes[i] is None:
      serr = None
    else:
      serr = subprocess.PIPE
    p = subprocess.Popen(c, cwd=localTempDir, stdin=sin,
                         stdout=sout, stderr=serr)
    if inPipes[i] is None:
      sin = None
    else:
      if not os.path.exists(inPipes[i]):
        raise IOError('Unable to locate inPipe file: %s for command %s'
                      % (inPipes[i], ' '.join(c)))
      sin = open(inPipes[i], 'r').read()
    if outPipes[i] is None:
      pout, perr = p.communicate(sin)
      handleReturnCode(p.returncode, cmds[i])
    else:
      with open(outPipes[i], 'w') as f:
        f.write(p.communicate(sin)[0])
      handleReturnCode(p.returncode, cmds[i])


def handleReturnCode(retcode, cmd):
  """ handle the return codes from runCommands
  """
  if not isinstance(retcode, int):
    raise TypeError('handleReturnCode takes an integer for '
                    'retcode, not a %s.' % retcode.__class__)
  if not isinstance(cmd, list):
    raise TypeError('handleReturnCode takes a list for '
                    'cmd, not a %s.' % cmd.__class__)
  if retcode:
    if retcode < 0:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s SIGNAL:%d' %(' '.join(cmd), -retcode))
    else:
      raise RuntimeError('Experienced an error while trying to execute: '
                         '%s retcode:%d' %(' '.join(cmd), retcode))


def createSequenceFile(sequences, tmpDir, filename='seq.fa'):
  """ given a dict (key is name, value is sequence) return path to temp file.
  """
  seqfile = os.path.join(tmpDir, filename)
  with open(seqfile, 'w') as f:
    for name in sequences:
      f.write('>%s\n%s' % (name, sequences[name]))
  return seqfile


def createAlignmentFile(alignments, tmpDir):
  """ given a list of alignments, return path to a temp file.
  """
  alnfile = os.path.join(tmpDir, 'aln.psl')
  with open(alnfile, 'w') as f:
    for a in alignments:
      if isinstance(a, str):
        f.write('%s\n' % a)
      elif isinstance(a, lib_filter.PslRow):
        f.write('%s\n' % a.pslString())
  return alnfile


def createBedFile(bedLines, name, tmpDir):
  """ given a list of string bed lines, return path to a temp file.
  """
  bedFile = os.path.join(tmpDir, name)
  with open(bedFile, 'w') as f:
    for b in bedLines:
        f.write('%s\n' % b)
  return bedFile


def bedLine(chrom, chromStart, chromEnd, name, score=None, strand=None,
            thickStart=None, thickEnd=None, itemRgb=None, blockCount=None,
            blockSizes=None, blockStarts=None):
  """ Give the fields, create a bed line string
  """

  s = ('%s %d %d %s'
       % (chrom, chromStart, chromEnd, name))
  if score is not None:
    for v in [strand, thickStart, thickEnd, itemRgb,
              blockCount, blockSizes, blockStarts]:
      assert(v is not None)
    s += (' %d %s %d %d %s %d %s %s'
          % (score, strand, thickStart, thickEnd, itemRgb, blockCount,
          blockSizes, blockStarts))
  return s


def simplePsl(strand, qSize, qStart, qEnd, tSize, tStart, tEnd,
              blockSizes, qStarts, tStarts, qName='query', tName='target'):
  """ Given a few of the fields, create a PslRow object.
  """
  line = ('%d %d %d %d %d %d %d %d %s %s %d %d %d %s %d %d %d %d %s %s %s'
          % (1, 0, 0, 0, 0, 0, 0, 0, strand, qName, qSize, qStart, qEnd,
             tName, tSize, tStart, tEnd, len(blockSizes),
             ','.join([str(b) for b in blockSizes]),
             ','.join([str(b) for b in qStarts]),
             ','.join([str(b) for b in tStarts]),
             ))
  return lib_filter.PslRow(line)

def numberOfUniqueTranscripts(transcripts):
  """ Given a list of transcripts, return the number of unique names.
  """
  names = set()
  for t in transcripts:
    if t not in names:
      names.add(t)
  return len(names)


def transcriptIsNonsense(t):
  """ check to see if the annotations contain a nonsense label
  """
  for a in t.annotations:
    if 'nonsense' in a.labels:
      return True
  return False


def transcriptHasOutOfFrame(t):
  """ check to see if the annotations contain an outOfFrame label
  """
  for a in t.annotations:
    if 'outOfFrame' in a.labels:
      return True
  return False


def outOfFrameCodonsThis(t):
  """ report the number of outOfFrameCodonsThis
  """
  for a in t.annotations:
    if 'outOfFrame' in a.labels:
      for l in a.labels:
        if l.startswith('outOfFrameCodonsThis'):
          return int(l.split('_')[1])
  return 0


def outOfFrameCodonsThem(t):
  """ report the number of outOfFrameCodonsThem
  """
  for a in t.annotations:
    if 'outOfFrame' in a.labels:
      for l in a.labels:
        if l.startswith('outOfFrameCodonsThem'):
          return int(l.split('_')[1])
  return 0


def frameShiftingIndels(t):
  """ report the number of frameShiftingIndels
  """
  for a in t.annotations:
    if 'outOfFrame' in a.labels:
      for l in a.labels:
        if l.startswith('frameShiftingIndels'):
          return int(l.split('_')[1])
  return 0


def transcriptHasMutations(t):
  """ check to see if the annotations contain mutation labels
  """
  for a in t.annotations:
    if 'nonsynon' in a.labels:
      return True
    if 'synon' in a.labels:
      return True
  return False


def nonSynonCount(t):
  """ count the number of nonsynon labels in the transcript annotations.
  """
  count = 0
  for a in t.annotations:
    for l in a.labels:
      if l == 'nonsynon':
        count += 1
  return count


def synonCount(t):
  """ count the number of synon labels in the transcript annotations.
  """
  count = 0
  for a in t.annotations:
    for l in a.labels:
      if l == 'synon':
        count += 1
  return count


def firstNonSynon(t):
  """ return the of the first nonsynon label encountered
  """
  for a in t.annotations:
    for i, l in enumerate(a.labels):
      if l == 'nonsynon':
        return a.labels[i + 1]
  return ''


def firstSynon(t):
  """ return the of the first synon label encountered
  """
  for a in t.annotations:
    for i, l in enumerate(a.labels):
      if l == 'synon':
        return a.labels[i + 1]
  return ''


class sequenceGetterTests(unittest.TestCase):
  def test_getSequences(self):
    """ getSequences must read a fasta and return a dict of Sequence objects.
    """
    sequences = {'chrA':
                   ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n'
                    'ACGT')
                  ,
                 'bannana.chr0':
                   ('ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACTTT\n'
                    'ACGTACGTACGTACGTACGTACGTACGTACG\n')
                  }
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('getSequences'))
    testFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(testFile)
    for name in sequences:
      seq = sequences[name]
      seq = seq.replace('\n', '').strip()
      self.assertTrue(name in seqDict)
      self.assertEqual(seqDict[name].getLength(), len(seq))
      self.assertEqual(seqDict[name].getSequence(), seq)
    for name in seqDict:
      self.assertTrue(name in sequences)
    self.addCleanup(removeDir, tmpDir)


class alignmentGetterTests(unittest.TestCase):
  def test_getAlignment(self):
    """ getAlignment must read a psl file and return a list of PslRow objects.
    """
    alignments = [
      '141 0 0 0 0 0 0 0 - ENSMUST00000178550.1 141 0 141 scaffold-11326 33702 21065 21206 1 141, 0, 21065,',
      '309 0 0 0 0 0 0 0 - ENSMUST00000179623.1 309 0 309 scaffold-1475 11716 9284 9593 1 309, 0, 9284,',
      '700 0 0 0 3 5 2 4 - ENSMUST00000179112.1 705 0 705 scaffold-189833 540197 335509 336213 4 14,129,140,417, 0,15,146,288, 335509,335523,335654,335796,',
                  ]
    fields = ['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert',
              'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName',
              'qSize', 'qStart', 'qEnd', 'tName', 'tSize',
              'tStart', 'tEnd', 'blockCount']
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('getAlignments'))
    testFile = createAlignmentFile(alignments, tmpDir)
    libAlignments = lib_filter.getAlignments(testFile)
    for i, a in enumerate(alignments, 0):
      data = a.split()
      for j, field in enumerate(fields, 0):
        self.assertEqual(data[j], str(getattr(libAlignments[i], field)))
      for j, field in enumerate(['blockSizes', 'qStarts', 'tStarts'], 18):
        self.assertEqual(
          data[j], ','.join(map(str, getattr(libAlignments[i], field))) + ',')
    self.addCleanup(removeDir, tmpDir)


class transcriptIteratorTests(unittest.TestCase):
  def test_transcriptIterator(self):
    """ tests the transcriptIterator function
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))

    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    names = ['ENSMUST00000178026.1', 'ENSMUST00000095795.4',
             'ENSMUST00000034053.5', 'ENSMUST00000034053.5',
             'ENSMUST00000034053.5']
    self.assertEqual(len(transcripts), 5)
    for attr, values in [('name', names),
                         ('thickStart', [2812370, 2812370, 466, 4903, 14291]),
                         ('thickEnd', [3038729, 3038729, 4248, 4996, 24866]),
                         ('itemRgb', ['128,0,0', '128,0,0',
                                      '128,0,0', '128,0,0','128,0,0',]),
                         ('chromosomeInterval',
                          [
          lib_filter.ChromosomeInterval('1', 2812346, 3113743, False),
          lib_filter.ChromosomeInterval('1', 2812346, 3113783, True),
          lib_filter.ChromosomeInterval('scaffold-100021', 466, 4248, False),
          lib_filter.ChromosomeInterval('scaffold-138877', 4903, 5091, False),
          lib_filter.ChromosomeInterval('scaffold-2051', 13759, 24866, False),
          ]
                          ),
                         ]:
      for i, n in enumerate(values, 0):
        self.assertEquals(getattr(transcripts[i], attr), n)
    # exons
    self.assertEqual(len(transcripts[0].exons), 9)
    self.assertEqual(transcripts[0].exons[0].chromosome, '1')
    self.assertEqual(transcripts[0].exons[0].start, 2812346)
    self.assertEqual(transcripts[0].exons[0].stop, 2812346 + 54)
    self.assertEqual(transcripts[0].exons[0].strand, False)
    self.assertEqual(transcripts[1].exons[0].chromosome, '1')
    self.assertEqual(transcripts[1].exons[0].start, 2812346)
    self.assertEqual(transcripts[1].exons[0].stop, 2812346 + 54)
    self.assertEqual(transcripts[1].exons[0].strand, True)
    # annotations
    testAnnots = [
      [],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('1', 2812346, 2812349, True),
        'ENSMUST00000095795.4', ['noStop']),
       lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('1', 3113780, 3113783, True),
        'ENSMUST00000095795.4', ['noStart']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-100021', 466, 469, False),
        'ENSMUST00000034053.5', ['noStop']),
       lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-100021', 4245, 4248, False),
        'ENSMUST00000034053.5', ['noStart']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-138877', 4903, 4906, False),
        'ENSMUST00000034053.5', ['noStop']),
       ],
      [lib_filter.TranscriptAnnotation(
        lib_filter.ChromosomeInterval('scaffold-2051', 24863, 24866, False),
        'ENSMUST00000034053.5', ['noStart']),
       ],
      ]
    self.assertEqual(len(transcripts), len(testAnnots))
    for i in xrange(0, len(transcripts)):
      try:
        self.assertEqual(transcripts[i].annotations, testAnnots[i])
      except AssertionError:
        for j in xrange(0, len(transcripts[i].annotations)):
          tra, tea = transcripts[i].annotations[j], testAnnots[i][j]
          print('  %d  [%s %s %d %d %s %s %s %s]\n  %d  [%s %s %d %d %s %s %s %s]'
                % (j, tra.name,
                   tra.chromosomeInterval.chromosome,
                   tra.chromosomeInterval.start,
                   tra.chromosomeInterval.stop,
                   tra.chromosomeInterval.strand,
                   tra.name, str(tra.labels), tra._itemRgb,
                   j, tea.name,
                   tea.chromosomeInterval.chromosome,
                   tea.chromosomeInterval.start,
                   tea.chromosomeInterval.stop,
                   tea.chromosomeInterval.strand,
                   tea.name, str(tra.labels), tea._itemRgb,
                   ))
        raise
    # Check print functions
    self.assertEquals(transcripts[0].bedString().split(), transcriptBedLines[0].split())
    self.assertEquals(transcripts[1].bedString().split(), transcriptBedLines[1].split())
    self.assertEquals(transcripts[1].annotations[0].bedString().split(), transcriptDetailsBedLines[0].split())
    self.assertEquals(transcripts[1].annotations[1].bedString().split(), transcriptDetailsBedLines[1].split())
    # Check sort function for transcripts
    transcripts.reverse()
    names.reverse()
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].name, names[i])
    transcripts.sort()
    names.reverse()  # ASSUME THAT NAMES AND BED LINES IN TEST ARE ALREADY IN SORTED ORDER
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].name, names[i])

  def test_transcriptWriter(self):
    """ transcripts should be written out correctly.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptWriter'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    # test equality.
    self.assertEquals(len(transcripts), len(writtenTranscripts))
    transcripts.sort()
    writtenTranscripts.sort()
    for i in xrange(0, len(transcripts)):
      self.assertEquals(transcripts[i].chromosomeInterval,
                        writtenTranscripts[i].chromosomeInterval)
      self.assertEquals(transcripts[i].annotations,
                        writtenTranscripts[i].annotations)
      self.assertEquals(transcripts[i], writtenTranscripts[i])
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_transcriptObjectManagement(self):
    """ Transcript and TranscriptAnnotation counts should not change.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113743, 'ENSMUST00000178026.1', 0, '-', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,90,165,105,13,45',
        '0,58,62,698,1209,1305,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        '1', 2812346, 3113783, 'ENSMUST00000095795.4', 0, '+', 2812370, 3038729,
        '128,0,0', 9, '54,2,89,249,197,52,105,13,85',
        '0,58,62,698,1209,1418,226292,301050,301352'))
    transcriptBedLines.append(bedLine(
        'scaffold-100021', 466, 4248, 'ENSMUST00000034053.5', 0, '-', 466, 4248,
        '128,0,0', 2, '85,152', '0,3630'))
    transcriptBedLines.append(bedLine(
        'scaffold-138877', 4903, 5091, 'ENSMUST00000034053.5', 0, '-', 4903,
        4996, '128,0,0', 1, '188', '0'))
    transcriptBedLines.append(bedLine(
        'scaffold-2051', 13759, 24866, 'ENSMUST00000034053.5', 0, '-', 14291,
        24866, '128,0,0', 4, '722,112,131,188', '0,1977,4316,10919'))
    # the following bed lines can cause problems because they represent
    # multilpe instances of a transcript on a single chromosome.
    transcriptBedLines += [
      'X 73726551 74227830 ENSMUST00000101560.3 0 + 73726551 74227830 128,0,0 12 110,22,6,14,4,15,21,18,143,6,21,23 0,500953,500979,500987,501006,501011,501028,501050,501069,501213,501220,501256',
      'X 73726990 74230505 ENSMUST00000101560.3 0 + 73726990 73727560 128,0,0 30 15,10,109,50,54,225,4,44,25,19,1,18,2,56,4,13,20,18,23,199,320,369,88,137,117,67,107,102,348,5 0,19,34,144,195,250,480,487,532,559,580,582,601,501428,501486,501491,501507,501529,501550,501576,501789,502110,502481,502571,502709,502827,502924,503055,503158,503510',
      'X 74211514 74218259 ENSMUST00000101560.3 0 + 74213333 74218259 128,0,0 8 40,15,18,8,148,254,70,133 0,1726,1742,1761,1770,1919,6540,6612',
      ]
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        '1', 2812346, 2812349, 'noStop/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        '1', 3113780, 3113783, 'noStart/ENSMUST00000095795.4'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 466, 469, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-100021', 4245, 4248, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-138877', 4903, 4906, 'noStop/ENSMUST00000034053.5'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-2051', 24863, 24866, 'noStart/ENSMUST00000034053.5'))
    transcriptDetailsBedLines += [
      'X 73726551 73726554 noStart/ENSMUST00000101560.3',
      'X 73726551 73726661 frameMismatch/ENSMUST00000101560.3',
      'X 73726551 74227830 badFrame/ENSMUST00000101560.3',
      'X 73726990 73726993 noStart/ENSMUST00000101560.3',
      'X 73726990 73727560 badFrame/ENSMUST00000101560.3',
      'X 73727005 73727009 cdsGap/ENSMUST00000101560.3',
      'X 73727009 73727019 frameDiscontig/ENSMUST00000101560.3',
      'X 73727019 73727024 cdsGap/ENSMUST00000101560.3',
      'X 73727133 73727134 cdsGap/ENSMUST00000101560.3',
      'X 73727134 73727184 frameDiscontig/ENSMUST00000101560.3',
      'X 73727184 73727185 cdsGap/ENSMUST00000101560.3',
      'X 73727185 73727239 frameDiscontig/ENSMUST00000101560.3',
      'X 73727239 73727240 cdsGap/ENSMUST00000101560.3',
      'X 73727465 73727470 cdsGap/ENSMUST00000101560.3',
      'X 73727470 73727474 frameDiscontig/ENSMUST00000101560.3',
      'X 73727474 73727477 cdsGap/ENSMUST00000101560.3',
      'X 73727521 73727522 cdsGap/ENSMUST00000101560.3',
      'X 73727522 73727547 frameDiscontig/ENSMUST00000101560.3',
      'X 73727547 73727549 cdsGap/ENSMUST00000101560.3',
      'X 73727559 73727560 noStop/ENSMUST00000101560.3',
      'X 73727568 73727570 utrGap/ENSMUST00000101560.3',
      'X 73727571 73727572 utrGap/ENSMUST00000101560.3',
      'X 73727590 73727591 utrGap/ENSMUST00000101560.3',
      'X 73727593 74228418 unknownUtrSplice/GA..AA/ENSMUST00000101560.3',
      'X 74213255 74213256 utrGap/ENSMUST00000101560.3',
      'X 74213274 74213275 utrGap/ENSMUST00000101560.3',
      'X 74213283 74213284 utrGap/ENSMUST00000101560.3',
      'X 74213333 74218259 badFrame/ENSMUST00000101560.3',
      'X 74213432 74213433 cdsGap/ENSMUST00000101560.3',
      'X 74213433 74213687 frameDiscontig/ENSMUST00000101560.3',
      'X 74218124 74218126 cdsGap/ENSMUST00000101560.3',
      'X 74218258 74218259 noStop/ENSMUST00000101560.3',
      'X 74227504 74227526 frameDiscontig/ENSMUST00000101560.3',
      'X 74227526 74227530 cdsGap/ENSMUST00000101560.3',
      'X 74227530 74227536 frameDiscontig/ENSMUST00000101560.3',
      'X 74227536 74227538 cdsGap/ENSMUST00000101560.3',
      'X 74227552 74227557 cdsGap/ENSMUST00000101560.3',
      'X 74227561 74227562 cdsGap/ENSMUST00000101560.3',
      'X 74227562 74227577 frameDiscontig/ENSMUST00000101560.3',
      'X 74227577 74227579 cdsGap/ENSMUST00000101560.3',
      'X 74227579 74227600 frameDiscontig/ENSMUST00000101560.3',
      'X 74227600 74227601 cdsGap/ENSMUST00000101560.3',
      'X 74227601 74227619 frameDiscontig/ENSMUST00000101560.3',
      'X 74227619 74227620 cdsGap/ENSMUST00000101560.3',
      'X 74227620 74227763 frameDiscontig/ENSMUST00000101560.3',
      'X 74227763 74227764 cdsGap/ENSMUST00000101560.3',
      'X 74227764 74227770 frameDiscontig/ENSMUST00000101560.3',
      'X 74227770 74227771 cdsGap/ENSMUST00000101560.3',
      'X 74227792 74227807 cdsGap/ENSMUST00000101560.3',
      'X 74227829 74227830 noStop/ENSMUST00000101560.3',
      'X 74228474 74228476 utrGap/ENSMUST00000101560.3',
      'X 74228480 74228481 utrGap/ENSMUST00000101560.3',
      'X 74228494 74228497 utrGap/ENSMUST00000101560.3',
      'X 74228517 74228519 utrGap/ENSMUST00000101560.3',
      'X 74228537 74228540 utrGap/ENSMUST00000101560.3',
      'X 74228563 74228566 utrGap/ENSMUST00000101560.3',
      'X 74228765 74228779 utrGap/ENSMUST00000101560.3',
      'X 74229099 74229100 utrGap/ENSMUST00000101560.3',
      'X 74229469 74229471 utrGap/ENSMUST00000101560.3',
      'X 74229559 74229561 utrGap/ENSMUST00000101560.3',
      'X 74229698 74229699 utrGap/ENSMUST00000101560.3',
      'X 74229816 74229817 utrGap/ENSMUST00000101560.3',
      'X 74229884 74229914 unknownUtrSplice/AA..GT/ENSMUST00000101560.3',
      'X 74230021 74230045 unknownUtrSplice/AA..GA/ENSMUST00000101560.3',
      'X 74230147 74230148 utrGap/ENSMUST00000101560.3',
      'X 74230496 74230500 utrGap/ENSMUST00000101560.3',
      ]
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    originalTranscriptCount = len(transcripts)
    originalTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), transcripts))
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptObjectManagement'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    newTranscriptCount = len(writtenTranscripts)
    newTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), writtenTranscripts))
    # test equality.
    self.assertEquals(originalTranscriptCount, newTranscriptCount)
    self.assertEquals(originalTranscriptAnnotationCount,
                      newTranscriptAnnotationCount)
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_transcriptReading_0(self):
    """ Transcript objects should be read in correctly
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'scaffold-444', 41415, 87033, 'ENSMUST00000169901.2', 0, '-', 41415,
        45156, '128,0,0', 5, '128,12,219,90,27', '0,131,3590,42232,45591'))
    transcriptBedLines.append(bedLine(
        'scaffold-444', 72633, 82553, 'ENSMUST00000169901.2', 0, '-', 72782,
        82485, '0,128,0', 5, '51,156,104,140,219', '0,129,4370,7482,9701'))
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41415, 41418,
        'noStop/alignmentPartialMap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41543, 41546,
        'cdsGap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 72633, 82553,
        'hasBadCopies/count_1/ENSMUST00000169901.2'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    originalTranscriptCount = len(transcripts)
    originalTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), transcripts))
    # test transcript count
    self.assertEquals(originalTranscriptCount, 2)
    self.assertEquals(len(transcripts[0].annotations), 3)
    self.assertEquals(len(transcripts[1].annotations), 1)
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('transcriptReading'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(outBed, outDetailsBed)
    newTranscriptCount = len(writtenTranscripts)
    newTranscriptAnnotationCount = sum(
      map(lambda t:len(t.annotations), writtenTranscripts))
    # test transcript count
    self.assertEquals(newTranscriptCount, 2)
    self.assertEquals(len(writtenTranscripts[0].annotations), 3)
    self.assertEquals(len(writtenTranscripts[1].annotations), 1)
    # test equality.
    self.assertEquals(originalTranscriptCount, newTranscriptCount)
    self.assertEquals(originalTranscriptAnnotationCount,
                      newTranscriptAnnotationCount)
    # cleanup
    self.addCleanup(removeDir, tmpDir)


class pslCoordinateSpaceTests(unittest.TestCase):
  def test_psl_targetCoordinateToQuery(self):
    """ PSLRow.targetCoordinateToQuery() should return correct information.
    """
    psls = []
    psls.append(simplePsl('+', 10, 0, 10, 10, 0, 10, [10], [0], [0]))
    psls.append(simplePsl('+', 10, 3, 10, 10, 3, 10, [7], [3], [3]))
    psls.append(simplePsl('+', 10, 3, 8, 10, 3, 8, [5], [3], [3]))
    psls.append(simplePsl('+', 10, 3, 8, 10, 1, 6, [5], [3], [1]))
    psls.append(simplePsl('+', 20, 3, 17, 30, 1, 22,
                          [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    psls.append(simplePsl('-', 24, 3, 17, 30, 1, 22,
                          [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    psls.append(simplePsl('+', 1257, 0, 683, 13608, 758, 4954,
                          [109, 90, 484], [0, 109, 199], [758,2563,4470]))
    psls.append(simplePsl('-', 5353, 101, 3546, 31353, 11338, 31141,
                          [168,136,213,816,100,571],
                          [1807,2922,3058,3271,4581,4681],
                          [11338,13222,13580,13796,29960,30570],))
    psls.append(simplePsl('+', 61, 22, 61, 61, 0, 61, [14, 5], [22, 56], [0, 14]))
    psls.append(simplePsl('-', 61, 4, 56, 61, 0, 61, [20, 18], [5, 39], [0, 20]))
    ####################
    #   0        9
    # q ++++++++++
    # t ++++++++++
    #   0        9
    psl = psls[0]
    for i in xrange(0, 10):
      self.assertEqual(i, psl.targetCoordinateToQuery(i))
    for i in xrange(10, 15):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    ####################
    #      3     9
    # q    +++++++
    # t    +++++++
    #      3     9
    psl = psls[1]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(3, 10):
      self.assertEqual(i, psl.targetCoordinateToQuery(i))
    for i in xrange(10, 15):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    ####################
    #      3   7
    # q    +++++
    # t    +++++
    #      3   7
    psl = psls[2]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(3, 8):
      self.assertEqual(i, psl.targetCoordinateToQuery(i))
    for i in xrange(8, 10):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    ####################
    #      3   7
    # q    +++++
    # t    +++++
    #      1   5
    psl = psls[3]
    for i in xrange(0, 1):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(1, 6):
      self.assertEqual(i + 2, psl.targetCoordinateToQuery(i))
    for i in xrange(6, 10):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    ####################
    #      3   7  8 10 15
    # q    +++++  +++  ++
    # t    +++++  +++  ++
    #      1   5  10   20
    psl = psls[4]
    for i in xrange(0, 1):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(1, 6):
      self.assertEqual(i + 2, psl.targetCoordinateToQuery(i))
    for i in xrange(6, 10):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(10, 13):
      self.assertEqual(i - 2, psl.targetCoordinateToQuery(i))
    for i in xrange(13, 20):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(20, 22):
      self.assertEqual(i - 5, psl.targetCoordinateToQuery(i))
    for i in xrange(22, 25):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    ####################
    # psls.append(simplePsl('-', 24, 3, 17, 30, 1, 22,
    #                      [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    #      1   5  10,12  20,21
    # t    +++++  +++    ++
    # q    -----  ---    -- (24 long)
    #     -3   7  8 10   15,16
    #     +20  16 15,13  8,7
    psl = psls[5]
    for q, t in [(None, 0), (20, 1), (16, 5), (None, 6),
                 (None, 7), (15, 10), (13, 12), (None, 13),
                 (None, 19), (8, 20), (7, 21), (None, 22)
                 ]:
      self.assertEqual(q, psl.targetCoordinateToQuery(t))
    ####################
    # psls.append(simplePsl('+', 1257, 0, 683,
    #                       13608, 758, 4954,
    #                       [109, 90, 484], [0, 109, 199], [758,2563,4470]))
    #
    #    + 0    108       109   198     199   682
    # q    ++++++         +++++++       +++++++
    # t    ++++++         +++++++       +++++++
    #    + 758  866       2563  2652    4470  4953
    psl = psls[6]
    for i in xrange(750, 758):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(4954, 4960):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for q, t in [(0, 758), (108, 866), (None, 867),
                 (None, 2562), (109, 2563), (198, 2652), (None, 2653),
                 (None, 4469), (199, 4470), (682, 4953), (None, 4954),
                 (130, 2584)
                 ]:
      self.assertEqual(q, psl.targetCoordinateToQuery(t))
    ####################
    # psls.append(simplePsl('-', 5353, 101, 3546, 31353, 11338, 31141,
    #                       [168,136,213,816,100,571],
    #                       [1807,2922,3058,3271,4581,4681],
    #                       [11338,13222,13580,13796,29960,30570],))
    #
    #
    #
    #  -
    #  +1807  1974    2922  3057    3058  3270    3271  4086    4581  4680    4681  5251
    # q -------       -------       -------       -------       -------       -------
    # t +++++++       +++++++       +++++++       +++++++       +++++++       +++++++
    #   11338 11505   13222 13557   13580 13792   13796 14611   29960 30059   30570 31140
    psl = psls[7]
    for i in xrange(11335, 11338):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for i in xrange(31141, 31145):
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for q, t in [(1807, 11338), (1974, 11505),
                 (None, 11506), (None, 13221),
                 (2922, 13222), (3057, 13557),
                 (None, 13558), (None, 13579),
                 (3058, 13580), (3270, 13792),
                 (None, 13793), (None, 13795),
                 (3171, 13796), (4086, 14611),
                 (None, 14612), (None, 29961),
                 (4581, 29960), (4680, 30059),
                 (None, 30060), (None, 30569),
                 (4681, 30570), (5251, 31140),]:
      pass
      #self.assertEqual(q, psl.targetCoordinateToQuery(t))
    ####################
    #0         1         2         3         4         5         6 tens position in query
    #0123456789012345678901234567890123456789012345678901234567890 ones position in query
    #                      ++++++++++++++                    +++++ plus strand alignment on query
    #    ------------------              --------------------      minus strand alignment on query
    #0987654321098765432109876543210987654321098765432109876543210 ones pos in query negative strand coordinates
    #6         5         4         3         2         1         0 tens pos in query negative strand coordinates
    # Plus strand:
    #      qStart=22
    #      qEnd=61
    #      blockSizes=14,5
    #      qStarts=22,56
    # Minus strand:
    #      qStart=4
    #      qEnd=56
    #      blockSizes=20,18
    #      qStarts=5,39
    # psls.append(simplePsl('+', 61, 22, 61, 61, 0, 61, [14, 5], [22, 56], [0, 14]))
    # psls.append(simplePsl('-', 61, 4, 56, 61, 0, 61, [20, 18], [5, 39], [0, 20]))
    psl = psls[8]
    for i in [19, 20, 21]:
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for q, t in [(22, 0), (35, 13),
                 (56, 14), (60, 18),
                 ]:
      self.assertEqual(q, psl.targetCoordinateToQuery(t))
    #     0   19  20   37
    # t   +++++   ++++++
    # q   -----   ------
    #    -5   24  39   56
    #    +55  36  21   4
    psl = psls[9]
    for i in [38, 39, 40]:
      self.assertEqual(None, psl.targetCoordinateToQuery(i))
    for q, t in [(55, 0), (36, 19),
                 (21, 20), (4, 37),
                 ]:
      self.assertEqual(q, psl.targetCoordinateToQuery(t))


  def test_psl_queryCoordinateToTarget(self):
    """ PSLRow.queryCoordinateToTarget() should return correct information.
    """
    psls = []
    psls.append(simplePsl('+', 10, 0, 10, 10, 0, 10, [10], [0], [0]))
    psls.append(simplePsl('+', 10, 3, 10, 10, 3, 10, [7], [3], [3]))
    psls.append(simplePsl('+', 10, 3, 8, 10, 3, 8, [5], [3], [3]))
    psls.append(simplePsl('+', 10, 3, 8, 10, 1, 6, [5], [3], [1]))
    psls.append(simplePsl('+', 20, 3, 17, 30, 1, 22,
                          [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    psls.append(simplePsl('-', 24, 3, 17, 30, 1, 22,
                          [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    ####################
    #   0        9
    # q ++++++++++
    # t ++++++++++
    #   0        9
    psl = psls[0]
    for i in xrange(0, 10):
      self.assertEqual(i, psl.queryCoordinateToTarget(i))
    for i in xrange(10, 15):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    ####################
    #      3     9
    # q    +++++++
    # t    +++++++
    #      3     9
    psl = psls[1]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    for i in xrange(3, 10):
      self.assertEqual(i, psl.queryCoordinateToTarget(i))
    for i in xrange(10, 15):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    ####################
    #      3   7
    # q    +++++
    # t    +++++
    #      3   7
    psl = psls[2]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    for i in xrange(3, 8):
      self.assertEqual(i, psl.queryCoordinateToTarget(i))
    for i in xrange(8, 10):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    ####################
    #      3   7
    # q    +++++
    # t    +++++
    #      1   5
    psl = psls[3]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    for i in xrange(3, 8):
      self.assertEqual(i - 2, psl.queryCoordinateToTarget(i))
    for i in xrange(8, 10):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    ####################
    #      3   7  8 10 15
    # q    +++++  +++  ++
    # t    +++++  +++  ++
    #      1   5  10   20
    psl = psls[4]
    for i in xrange(0, 3):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    for i in xrange(3, 8):
      self.assertEqual(i - 2, psl.queryCoordinateToTarget(i))
    for i in xrange(8, 11):
      self.assertEqual(i + 2, psl.queryCoordinateToTarget(i))
    for i in xrange(11, 15):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    for i in xrange(15, 17):
      self.assertEqual(i + 5, psl.queryCoordinateToTarget(i))
    for i in xrange(17, 20):
      self.assertEqual(None, psl.queryCoordinateToTarget(i))
    ####################
    #      1   5  10,12  20,21
    # t    +++++  +++    ++
    # q    -----  ---    -- (24 long)
    #     -3   7  8 10   15,16
    #     +20  16 15,13  8,7
    psl = psls[5]
    for q, t in [(21, None), (20, 1), (16, 5),
                 (15, 10), (13, 12), (12, None),
                 (9, None), (8, 20), (7, 21), (6, None),
                 ]:
      self.assertEqual(t, psl.queryCoordinateToTarget(q))
  def test_psl_roundtrip_queryTarget(self):
    """ PSLRow.queryCoordinateToTarget() and PSLRow.targetCoordinateToQuery() should play nice.
    """
    psls = []
    psls.append(simplePsl('-', 24, 3, 17, 30, 1, 22,
                          [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    # psls.append(simplePsl('-', 24, 3, 17, 30, 1, 22,
    #                      [5, 3, 2], [3, 8, 15], [1, 10, 20]))
    #      1   5  10,12  20,21
    # t    +++++  +++    ++
    # q    -----  ---    -- (24 long)
    #     -3   7  8 10   15,16
    #     +20  16 15,13  8,7
    psl = psls[0]
    for q, t in [(20, 1), (16, 5),
                 (15, 10), (13, 12),
                 (8, 20), (7, 21),
                 ]:
      self.assertEqual(t, psl.queryCoordinateToTarget(psl.targetCoordinateToQuery(t)))
      self.assertEqual(q, psl.targetCoordinateToQuery(psl.queryCoordinateToTarget(q)))


class codonGeneSpaceTests(unittest.TestCase):
  def test_transcript_getMRna_0(self):
    """ Transcript.getMRna() should return correct information.
    """
    seq = lib_filter.Sequence(
    #             0        9          20        30        40
      'test_a', 'NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn')
    #              *********              *********
    #              0        9          20        30
    seq_rc = lib_filter.Sequence(
    #             0        9          20        30        40
      'test_rc', 'nnnNNNNNNCTACTCCgCCTnnnnnnnnnnACGaGaaaCATNN')
    #                      *********              *********
    #                      0        9          20 23
    #                      GATGAGGcGGAnnnnnnnnnnTGCtCtttGTA
    seq_thickThin = lib_filter.Sequence(
      'test_thickThin',
    #  0        9          20        30        40        50        60        70
      'nnnACGTACGTACGTACGTAACTACGTACGttATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn\n')
    #     -----------------------------*********              *********
    #     tttttttttttttttttttttttttttttTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    #     0        9          20        30        40        50        60        70
    seq_thickThin_splitExon_0 = lib_filter.Sequence(
      'test_thickThin_splitExon_0',
    #  0        9          20        30        40        50        60        70
      'nnnACGTACGTACGTACGTAACTACGTACGttATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn\n')
    #     -----                        *********              *********
    #     ttttttttttttttttttttttttTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    #     0        9          20        30        40        50        60        70
    seq_thickThin_splitExon_1 = lib_filter.Sequence(
      'test_thickThin_splitExon_1',
    #  0        9          20        30        40        50        60        70
      'nnnACGTACGTACGTACGTAACTACGTACGttATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn\n')
    #     -----                   -----*********              *********-----
    #     tttttttttttttttttttttttttttttTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTttttt
    #     0        9          20        30        40        50        60        70
    seq_thickThin_splitExon_2 = lib_filter.Sequence(
      'test_thickThin_splitExon_2',
    #  0        9          20        30        40        50        60        70
      'nnnNNNNNNCTACTCCgCCTnnnnnnnnnnACGaGaaaCATaaCGTACGTAGTTACGTACGTACGTACGTnnn\n')
    #      -----*********              *********-----                   -----
    #      tttttTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTttttttttttttttttttttttttttttt
    #      0        9          20        30        40        50        60        70

    truth = 'ATGtttCtCGcGGAGTAG'
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_rc', 9, 41, 'gene', 0, '-', 9, 41,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_thickThin', 3, 64, 'gene', 0, '+', 32, 64,
        '128,0,0', 2, '38,9',
        '0,52'))
    transcriptBedLines.append(bedLine(
        'test_thickThin_splitExon_0', 3, 64, 'gene', 0, '+', 32, 64,
        '128,0,0', 3, '5,9,9',
        '0,29,52'))
    transcriptBedLines.append(bedLine(
        'test_thickThin_splitExon_1', 3, 69, 'gene', 0, '+', 32, 64,
        '128,0,0', 3, '5,14,14',
        '0,24,52'))
    transcriptBedLines.append(bedLine(
        'test_thickThin_splitExon_2', 4, 70, 'gene', 0, '-', 9, 41,
        '128,0,0', 3, '14,14,5',
        '0,28,61'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    mrna = transcripts[0].getMRna(seq)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)
    mrna = transcripts[1].getMRna(seq_rc)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)
    mrna = transcripts[2].getMRna(seq_thickThin)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)
    mrna = transcripts[3].getMRna(seq_thickThin_splitExon_0)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)
    mrna = transcripts[4].getMRna(seq_thickThin_splitExon_1)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)
    mrna = transcripts[5].getMRna(seq_thickThin_splitExon_2)
    self.assertEqual(len(truth), len(mrna))
    self.assertEqual(truth, mrna)

  def test_transcript_mRnaCoordinateToExon(self):
    """ mRnaCoordinateToExon() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    for i in xrange(0, 6):
      self.assertEqual(2 + i, transcripts[0].mRnaCoordinateToExon(i))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    for i in xrange(0, 6):
      self.assertEqual(2 + i, transcripts[1].mRnaCoordinateToExon(i))
    ##############################
    # real data, ENSMUST00000179734.1-0
    # C382543 0 2237 ENSMUST00000179734.1-0 0 + 618 2074 128,0,0 6 381,20,826,288,74,97 0,392,472,1317,2047,2140
    for i in xrange(0, 680):
      # 547 = 381 + 20 + 146
      self.assertEqual(547 + i, transcripts[2].mRnaCoordinateToExon(i))
    for i in xrange(0, 288):
      # 1227 = 381 + 20 + 146 + 680
      self.assertEqual(1227 + i, transcripts[2].mRnaCoordinateToExon(680 + i))
    for i in xrange(0, 27):
      # 1505 = 381 + 20 + 146 + 680 + 288
      self.assertEqual(1515 + i, transcripts[2].mRnaCoordinateToExon(680 + 288 + i))

  def test_transcript_mRnaCoordinateToChromosome(self):
    """ mRnaCoordinateToChromosome() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    for i in xrange(0, 2):
      self.assertEqual(3 + i, transcripts[0].mRnaCoordinateToChromosome(i))
    for i in xrange(2, 6):
      self.assertEqual(4 + i, transcripts[0].mRnaCoordinateToChromosome(i))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    for i in xrange(0, 2):
      self.assertEqual(9 - i, transcripts[1].mRnaCoordinateToChromosome(i))
    for i in xrange(2, 6):
      self.assertEqual(6 + 2 - i, transcripts[1].mRnaCoordinateToChromosome(i))
    ##############################
    # real data, ENSMUST00000179734.1-0
    # C382543 0 2237 ENSMUST00000179734.1-0 0 + 618 2074 128,0,0 6 381,20,826,288,74,97 0,392,472,1317,2047,2140
    for i in xrange(0, 680):
      # 618 is the max of thick start and exon 3
      self.assertEqual(618 + i, transcripts[2].mRnaCoordinateToChromosome(i))
    for i in xrange(0, 288):
      # 1317 = start of exon 4
      self.assertEqual(1317 + i, transcripts[2].mRnaCoordinateToChromosome(680 + i))
    for i in xrange(0, 27):
      # 2047 is start of exon 5
      self.assertEqual(2047 + i, transcripts[2].mRnaCoordinateToChromosome(680 + 288 + i))



  def test_transcript_mRnaCoordinateToCodon(self):
    """ mRnaCoordinateToCodon() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    values = [(0, (0, 0)),
              (1, (0, 1)),
              (2, (0, 2)),
              (3, (1, 0)),
              (4, (1, 1)),
              (5, (1, 2)),
              (6, (2, 0)),
              (12, (4, 0)),
              (13, (4, 1)),
              (14, (4, 2)),
              ]
    t = transcripts[0]
    for grain, flour in values:
      self.assertEqual(flour, t.mRnaCoordinateToCodon(grain))

  def test_transcript_codonCoordinateToMRna(self):
    """ codonCoordinateToMRna() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    values = [(0, (0, 0)),
              (1, (0, 1)),
              (2, (0, 2)),
              (3, (1, 0)),
              (4, (1, 1)),
              (5, (1, 2)),
              (6, (2, 0)),
              (12, (4, 0)),
              (13, (4, 1)),
              (14, (4, 2)),
              ]
    t = transcripts[0]
    for flour, grain in values:
      self.assertEqual(flour, t.codonCoordinateToMRna(grain))

  def test_transcript_roundTripCodonMRna(self):
    """ mRnaCordinateToCodon and codonCoordinateToMRna should play nice.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    values = [(0, (0, 0)),
              (1, (0, 1)),
              (2, (0, 2)),
              (3, (1, 0)),
              (4, (1, 1)),
              (5, (1, 2)),
              (6, (2, 0)),
              (12, (4, 0)),
              (13, (4, 1)),
              (14, (4, 2)),
              ]
    t = transcripts[0]
    for flour, grain in values:
      self.assertEqual(flour, t.codonCoordinateToMRna(grain))
      self.assertEqual(grain, t.mRnaCoordinateToCodon(flour))
      self.assertEqual(grain,
                       t.mRnaCoordinateToCodon(t.codonCoordinateToMRna(grain)))
      self.assertEqual(flour,
                       t.codonCoordinateToMRna(t.mRnaCoordinateToCodon(flour)))

  def test_transcript_exonCoordinateToMRna(self):
    """ exonCoordinateToMRna() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptBedLines.append(bedLine(
        '4', 42458750, 42459176, 'ENSMUST00000098118.1', 0, '-', 42458750,
        42459176, '0,128,0', 1, '426', '0'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    t = transcripts[0]
    for i in xrange(0, 2):
      self.assertEqual(None, t.exonCoordinateToMRna(i))
    for i in xrange(2, 8):
      self.assertEqual(i - 2, t.exonCoordinateToMRna(i))
    self.assertEqual(None, t.exonCoordinateToMRna(8))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    t = transcripts[1]
    for i in xrange(0, 2):
      self.assertEqual(None, t.exonCoordinateToMRna(i))
    for i in xrange(2, 8):
      self.assertEqual(i - 2, t.exonCoordinateToMRna(i))
    self.assertEqual(None, t.exonCoordinateToMRna(8))
    self.assertEqual(None, t.exonCoordinateToMRna(None))
    ##############################
    # 4 42458750  42459176  ENSMUST00000098118.1  0 - 42458750  42459176  0,128,0 1 426 0
    # negative strand
    #              [0     426)
    #               |     |
    # mrna          +++++++
    # exon          +++++++
    #               |     |
    #              [0     426)
    # chromosome nnnnnnnnnnnnn
    #               |     |
    #       (42459176     42458750]
    t = transcripts[4]
    for i in xrange(0, 426):
      self.assertEqual(i, t.exonCoordinateToMRna(i))
    for i in xrange(426, 430):
      self.assertEqual(None, t.exonCoordinateToMRna(i))
  def test_transcript_roundtripExonMRna(self):
    """ exonCoordinateToMRna() and mRnaCoordinateToExon() should play nice.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    t = transcripts[0]
    for i in xrange(0, 2):
      self.assertEqual(None, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(i)))
    for i in xrange(2, 8):
      self.assertEqual(i, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(i)))
    self.assertEqual(None, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(8)))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    t = transcripts[1]
    for i in xrange(0, 2):
      self.assertEqual(None, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(i)))
    for i in xrange(2, 8):
      self.assertEqual(i, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(i)))
    self.assertEqual(None, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(8)))
    self.assertEqual(None, t.mRnaCoordinateToExon(t.exonCoordinateToMRna(None)))

  def test_transcript_exonCoordinateToChromosome(self):
    """ exonCoordinateToChromosome() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_rc', 9, 40, 'gene', 0, '-', 9, 40,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # exon coordinates
    #   0       8              9       17
    #   |       |              |       |
    # NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn
    # |    |    |    |    |    |    |    |    | |
    # 0    5    10        20        30          42
    # chromosome coordinates
    for i in xrange(0, 9):
      self.assertEqual(2 + i, transcripts[0].exonCoordinateToChromosome(i))
    for i in xrange(9, 18):
      self.assertEqual(25 - 9 + i, transcripts[0].exonCoordinateToChromosome(i))
    ##############################
    # negative strand
    # exon coordinates
    #   0       8              9       17
    #   |       |              |       |
    # NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn
    # | |    |    |    |    |    |    |    |    |
    #   40        30        20        10   5    0
    # chromosome coordinates
    for i in xrange(0, 9):
      self.assertEqual(40 - i, transcripts[1].exonCoordinateToChromosome(i))
    for i in xrange(9, 18):
      self.assertEqual(17 + 9 - i, transcripts[1].exonCoordinateToChromosome(i))
    ##############################
    # real data, ENSMUST00000179734.1-0
    # C382543 0 2237 ENSMUST00000179734.1-0 0 + 618 2074 128,0,0 6 381,20,826,288,74,97 0,392,472,1317,2047,2140
    for i in xrange(0, 381):
      self.assertEqual(i, transcripts[2].exonCoordinateToChromosome(i))
    for i in xrange(0, 20):
      # 392 is start of exon 2
      self.assertEqual(392 + i, transcripts[2].exonCoordinateToChromosome(381 + i))
    for i in xrange(0, 826):
      # 472 is start of exon 2
      self.assertEqual(472 + i, transcripts[2].exonCoordinateToChromosome(381 + 20 + i))
    for i in xrange(0, 288):
      # 1317 is start of exon 3
      self.assertEqual(1317 + i, transcripts[2].exonCoordinateToChromosome(381 + 20 + 826 + i))
    for i in xrange(0, 74):
      # 2047 is start of exon 4
      self.assertEqual(2047 + i, transcripts[2].exonCoordinateToChromosome(381 + 20 + 826 + 288 + i))
    for i in xrange(0, 97):
      # 2140 is start of exon 4
      self.assertEqual(2140 + i, transcripts[2].exonCoordinateToChromosome(381 + 20 + 826 + 288 + 74 + i))

  def test_transcript_chromosomeCoordinateToExon(self):
    """ chromosomeCoordinateToExon() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    t = transcripts[0]
    self.assertEqual(None, t.chromosomeCoordinateToExon(0))
    for i in xrange(1, 5):
      self.assertEqual(i - 1, t.chromosomeCoordinateToExon(i))
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(5))
    for i in xrange(6, 11):
      self.assertEqual(i - 2, t.chromosomeCoordinateToExon(i))
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(11))
    self.assertEqual(None, t.chromosomeCoordinateToExon(12))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    t = transcripts[1]
    for i in xrange(0, 2):
      self.assertEqual(None, t.chromosomeCoordinateToExon(i))
    for i in xrange(2, 7):
      self.assertEqual(8 + 2 - i, t.chromosomeCoordinateToExon(i))
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(7))
    for i in xrange(8, 12):
      self.assertEqual(8 + 3 - i, t.chromosomeCoordinateToExon(i))
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(12))
  def test_transcript_roundtripChromosomeExon(self):
    """ chromosomeCoordinateToExon() and exonCoordinateToChromosome() should play nice.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    t = transcripts[0]
    self.assertEqual(None, t.chromosomeCoordinateToExon(0))
    for i in xrange(1, 5):
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
      self.assertEqual(i,
                       t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(5))
    for i in xrange(6, 9):
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
      self.assertEqual(t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i)),
                       t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i))
                       )
    for i in xrange(9, 11):
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    t = transcripts[1]
    for i in xrange(0, 2):
      self.assertEqual(None, t.chromosomeCoordinateToExon(i))
    for i in xrange(2, 7):
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
      self.assertEqual(i,
                       t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(7))
    for i in xrange(8, 12):
      self.assertEqual(i,
                       t.exonCoordinateToChromosome(t.chromosomeCoordinateToExon(i))
                       )
    self.assertEqual(None, t.chromosomeCoordinateToExon(12))
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_rc', 9, 40, 'gene', 0, '-', 9, 40,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # exon coordinates
    #   0       8              9       17
    #   |       |              |       |
    # NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn
    # |    |    |    |    |    |    |    |    | |
    # 0    5    10        20        30          42
    # chromosome coordinates
    t = transcripts[0]
    for i in xrange(0, 9):
      self.assertEqual(i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i)))
    for i in xrange(9, 18):
      self.assertEqual(i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i)))
    ##############################
    # negative strand
    # exon coordinates
    #   0       8              9       17
    #   |       |              |       |
    # NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn
    # | |    |    |    |    |    |    |    |    |
    #   40        30        20        10   5    0
    # chromosome coordinates
    t = transcripts[1]
    for i in xrange(0, 9):
      self.assertEqual(i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i)))
    for i in xrange(9, 18):
      self.assertEqual(i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i)))
    ##############################
    # real data, ENSMUST00000179734.1-0
    # C382543 0 2237 ENSMUST00000179734.1-0 0 + 618 2074 128,0,0 6 381,20,826,288,74,97 0,392,472,1317,2047,2140
    t = transcripts[2]
    for i in xrange(0, 381):
      self.assertEqual(i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(i)))
    for i in xrange(0, 20):
      # 392 is start of exon 2
      self.assertEqual(381 + i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(381 + i)))
    for i in xrange(0, 826):
      # 472 is start of exon 2
      self.assertEqual(381 + 20 + i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(381 + 20 + i)))
    for i in xrange(0, 288):
      # 1317 is start of exon 3
      self.assertEqual(381 + 20 + 826 + i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(381 + 20 + 826 + i)))
    for i in xrange(0, 74):
      # 2047 is start of exon 4
      self.assertEqual(381 + 20 + 826 + 288 + i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(381 + 20 + 826 + 288 + i)))
    for i in xrange(0, 97):
      # 2140 is start of exon 4
      self.assertEqual(381 + 20 + 826 + 288 + 74 + i, t.chromosomeCoordinateToExon(t.exonCoordinateToChromosome(381 + 20 + 826 + 288 + 74 + i)))

  def test_transcript_chromosomeCoordinateToMRna(self):
    """ chromosomeCoordinateToMRna() must return correct values.
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 1, 11, 'gene', 0, '+', 3, 10,
        '128,0,0', 2, '4,5',
        '0,5'))
    transcriptBedLines.append(bedLine(
        'test_rc', 2, 12, 'gene', 0, '-', 3, 10,
        '128,0,0', 2, '5,4',
        '0,6'))
    transcriptBedLines.append(bedLine(
        'C382543', 0, 2237, 'ENSMUST00000179734.1-0', 0, '+', 618, 2074,
        '128,0,0', 6, '381,20,826,288,74,97', '0,392,472,1317,2047,2140'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    ##############################
    # positive strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #            |    |    |
    #            0    5    10
    # so to go from mrna to exon, we must add on the difference
    # between the thick start and thin start from the "start".
    t = transcripts[0]
    for i in [0, 1, 2, 5, 10, 11, 12]:
      self.assertEqual(None, t.chromosomeCoordinateToMRna(i))
    for i in xrange(3, 5):
      self.assertEqual(i - 3, t.chromosomeCoordinateToMRna(i))
    for i in xrange(6, 10):
      self.assertEqual(i - 4, t.chromosomeCoordinateToMRna(i))
    ##############################
    # negative strand
    #               0     5
    #               |     |
    # mrna          ++ ++++
    # exon        ..++ ++++.  two exons (thick and thin parts)
    #             |     |  |
    #             0     5  8
    # chromosome nnnnnnnnnnnnn
    #              |    |    |
    #              10   5    0
    t = transcripts[1]
    for i in [0, 1, 2, 7, 10, 11, 12]:
      self.assertEqual(None, t.chromosomeCoordinateToMRna(i))
    for i in xrange(3, 7):
      self.assertEqual(5 + 3 - i, t.chromosomeCoordinateToMRna(i))
    for i in xrange(8, 10):
      self.assertEqual(8 + 1 - i, t.chromosomeCoordinateToMRna(i))

class codonAminoAcidTests(unittest.TestCase):
  def test_codonToAminoAcid(self):
    """ codonToAminoAcid() needs to return correct amino acids for all codons.
    """
    # these pairs were input separately from pairs in lib_filter,
    # (no copy-paste!)
    knownPairs = [('Ala', 'GCT,GCC,GCA,GCG,GCN'),
                  ('Arg', 'CGT,CGC,CGA,CGG,AGA,AGG,CGN,MGR'),
                  ('Asn', 'AAT,AAC,AAY'),
                  ('Asp', 'GAT,GAC,GAY'),
                  ('Cys', 'TGT,TGC,TGY'),
                  ('Gln', 'CAA,CAG,CAR'),
                  ('Glu', 'GAA,GAG,GAR'),
                  ('Gly', 'GGT,GGC,GGA,GGG,GGN'),
                  ('His', 'CAT,CAC,CAY'),
                  ('Ile', 'ATT,ATC,ATA,ATH'),
                  ('Leu', 'TTA,TTG,CTT,CTC,CTA,CTG,YTR,CTN'),
                  ('Lys', 'AAA,AAG,AAR'),
                  ('Met', 'ATG'),
                  ('Phe', 'TTT,TTC,TTY'),
                  ('Pro', 'CCT,CCC,CCA,CCG,CCN'),
                  ('Ser', 'TCT,TCC,TCA,TCG,AGT,AGC,TCN,AGY'),
                  ('Thr', 'ACT,ACC,ACA,ACG,ACN'),
                  ('Trp', 'TGG'),
                  ('Tyr', 'TAT,TAC,TAY'),
                  ('Val', 'GTT,GTC,GTA,GTG,GTN'),
                  ('Stop', 'TAA,TGA,TAG,TAR,TRA'),
                  ('???', 'LOL,OMG,IDK,WTF,,too long')
                  ]
    for aa, codons in knownPairs:
      for c in codons.split(','):
        self.assertEqual(aa, lib_filter.codonToAminoAcid(c))
        self.assertEqual(aa, lib_filter.codonToAminoAcid(c.lower()))
        if len(c) > 2:
          self.assertEqual(aa, lib_filter.codonToAminoAcid(c[0].lower() + c[1:]))
          self.assertEqual(aa, lib_filter.codonToAminoAcid(
              c[0] + c[1].lower() + c[2]))
          self.assertEqual(aa, lib_filter.codonToAminoAcid(c[0:2] + c[2].lower()))
          self.assertEqual(aa, lib_filter.codonToAminoAcid(c[0:2].lower() + c[2]))
          self.assertEqual(aa, lib_filter.codonToAminoAcid(c[0] + c[1:3].lower()))
          self.assertEqual(aa, lib_filter.codonToAminoAcid(
              c[0].lower() + c[1] + c[2].lower()))


class filterTests(unittest.TestCase):
  def test_uniquify_0(self):
    """ uniquify should make unique names for transcripts with identical names
    """
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'scaffold-444', 41415, 87033, 'ENSMUST00000169901.2', 0, '-', 41415,
        45156, '128,0,0', 5, '128,12,219,90,27', '0,131,3590,42232,45591'))
    transcriptBedLines.append(bedLine(
        'scaffold-444', 72633, 82553, 'ENSMUST00000169901.2', 0, '-', 72782,
        82485, '0,128,0', 5, '51,156,104,140,219', '0,129,4370,7482,9701'))
    transcriptBedLines.append(bedLine(
        'scaffold-banana', 72633, 82553, 'ENSMUST00000169901.2', 0, '-', 72782,
        82485, '0,128,0', 5, '51,156,104,140,219', '0,129,4370,7482,9701'))
    transcriptDetailsBedLines = []
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41415, 41418,
        'noStop/alignmentPartialMap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 41543, 41546,
        'cdsGap/hasOkCopies/count_1/ENSMUST00000169901.2'))
    transcriptDetailsBedLines.append(bedLine(
        'scaffold-444', 72633, 82553,
        'hasBadCopies/count_1/ENSMUST00000169901.2'))
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('uniquify'))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run uniquify
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    originalGeneCheckBed = os.path.join(tmpDir, 'dummy.txt')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'dummy.txt')
    alignment = os.path.join(tmpDir, 'dummy.txt')
    sequence = os.path.join(tmpDir, 'dummy.txt')
    referenceSequence = os.path.join(tmpDir, 'dummy.txt')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'uniquify', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))
    # test equality.
    self.assertEquals(numberOfUniqueTranscripts(writtenTranscripts), 3)
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_nonsense_0(self):
    """ nonsense should detect nonsense codons.
    """
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('nonsense'))
    sequences = {'test_a':  # has nonsense
                   'NNATGtttCtCGTnnnnnnnnnnAGtaaGAGTAGNNNNNNnnn\n',
                 'test_ok':  # is ok
                   'NNATGtttCtCGTnnnnnnnnnnAGGcGGAGTAGNNNNNNnnn\n',
                 'test_ok_rc':  # is ok
                   'nnnNNNNNNCTACTCcccCTnnnnnnnnnnACGaGaaaCATNN\n',
                 'test_ok_thickThin':  # is ok
                   ('nnnACGTACG'
                    'TACGTACGTA'
                    'ACTACGTACG'
                    'ttATGtttCt'
                    'CGTnnnnnnn'
                    'nnnAGGcGGA'
                    'GTAGNNNNNN'
                    'nnn\n'),
                 'test_rc':  # has nonsense
                   'nnnNNNNNNCTACTCctaCTnnnnnnnnnnACGaGaaaCATNN\n',
                  }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_a', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_ok', 2, 34, 'gene', 0, '+', 2, 34,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_ok_rc', 9, 41, 'gene', 0, '-', 9, 41,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptBedLines.append(bedLine(
        'test_ok_thickThin', 3, 64, 'gene', 0, '+', 32, 64,
        '128,0,0', 2, '39,9',
        '0,53'))
    transcriptBedLines.append(bedLine(
        'test_rc', 9, 41, 'gene', 0, '-', 9, 41,
        '128,0,0', 2, '9,9',
        '0,23'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run nonsense
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    originalGeneCheckBed = os.path.join(tmpDir, 'dummy.txt')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'dummy.txt')
    alignment = os.path.join(tmpDir, 'dummy.txt')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'dummy.txt')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'nonsense', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))
    # test equality.
    self.assertTrue(transcriptIsNonsense(writtenTranscripts[0]))
    self.assertFalse(transcriptIsNonsense(writtenTranscripts[1]))
    self.assertFalse(transcriptIsNonsense(writtenTranscripts[2]))
    self.assertFalse(transcriptIsNonsense(writtenTranscripts[3]))
    self.assertTrue(transcriptIsNonsense(writtenTranscripts[4]))
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_mRnaCompare_0(self):
    """ mRnaCompare should detect outOfFrame mRNAs.
    """
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('mRnaCompare_0'))
    sequences = {'test_0_nr':  # non-ref
                   'ATGATCCAATGA\n',  # 12
                 }
    refSequences = {'test_0_r':  # ref
                    'ATGACCTCCAAATGA\n',  # 15
                    'test_0_r_rc':
                      'TCATTTGGAGGTCAT\n',  # 15
                    }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0  3  6 8 10  14
    # ref        ATGACCTCCAAATGA  query
    #            |||         |||  in-frame codon alignments
    # non ref    ATGA--TCC-AATGA  target
    #            0  3  4 6 7   11
    #               =========     out of frame
    # number of out of frame codons wrt target: 2
    # number of out of frame codons wrt query: 3
    # number of frame shifting indels: 2
    #####
    #            14 11 8 6 4   0
    # ref        ATGACCTCCAAATGA  query
    #            |||         |||  in-frame codon alignments
    # non ref    ATGA--TCC-AATGA  target
    #            0  3  4 6 7   11
    #               =========     out of frame
    # number of out of frame codons wrt target: 2
    # number of out of frame codons wrt query: 3
    # number of frame shifting indels: 2
    pslLines = [simplePsl('+', 15, 0, 15, 12, 0, 12,
                          [4, 3, 5], [0, 6, 10], [0, 4, 7],
                          qName='ensmust0', tName='test_0_nr'),
                simplePsl('+', 15, 0, 15, 12, 0, 12,
                          [4, 3, 5], [0, 6, 10], [0, 4, 7],
                          qName='ensmust1', tName='test_0_nr')
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 0, 15, 'ensmust0', 0, '+', 0, 15,
        '128,0,0', 1, '15', '0'))
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 0, 15, 'ensmust1', 0, '-', 0, 15,
        '128,0,0', 1, '15', '0'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 0, 12, 'ensmust0', 0, '+', 0, 12,
        '128,0,0', 1, '12', '0'))
    transcriptBedLines.append(bedLine(
        'test_0_nr', 0, 12, 'ensmust1', 0, '-', 0, 12,
        '128,0,0', 1, '12', '0'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run mRnaCompare
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'mRnaCompare', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir,
      extra='--allowSingleExons')
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))
    # test equality.
    wt = writtenTranscripts[1]
    self.assertTrue(transcriptHasOutOfFrame(wt))
    self.assertEqual(2, outOfFrameCodonsThis(wt))
    self.assertEqual(3, outOfFrameCodonsThem(wt))
    self.assertEqual(2, frameShiftingIndels(wt))
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_mRnaCompare_1(self):
    """ mRnaCompare should detect outOfFrame mRNAs.
    """
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('mRnaCompare_1'))
    sequences = {'test_0_nr':  # non-ref / target
                   'ATGATTAAATGA\n',  # 12
                 }
    refSequences = {'test_0_r':  # ref / query
                    'ATGATCCAATGA\n'  # 12
                    }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0          11
    # ref        ATGATCCAATGA  query
    #            |||||**|||||  in-frame codon alignments
    # non ref    ATGATTAAATGA  target
    #            0          11
    #                  ***    nonsynon AAA.Lys_CAA.Gln
    #               ===       synon    ATT.Ile_ATC.Ile
    # number of out of frame codons wrt target: 0
    # number of frame shifting indels: 0
    #####
    pslLines = [simplePsl('+', 12, 0, 12, 12, 0, 12,
                          [12], [0], [0],
                          qName='ensmust0', tName='test_0_nr'),
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 0, 12, 'ensmust0', 0, '+', 0, 12,
        '128,0,0', 1, '12', '0'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 0, 12, 'ensmust0', 0, '+', 0, 12,
        '128,0,0', 1, '12', '0'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run mRnaCompare
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'mRnaCompare', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir,
      extra='--allowSingleExons')
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))
    # test equality.
    self.assertTrue(transcriptHasMutations(writtenTranscripts[0]))
    self.assertEqual(1, nonSynonCount(writtenTranscripts[0]))
    self.assertEqual(1, synonCount(writtenTranscripts[0]))
    self.assertEqual('AAA.Lys_CAA.Gln', firstNonSynon(writtenTranscripts[0]))
    self.assertEqual('ATT.Ile_ATC.Ile', firstSynon(writtenTranscripts[0]))
    # cleanup
    self.addCleanup(removeDir, tmpDir)

  def test_indel_0(self):
    """indel should produce deletion annotations correctly."""
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('indel_0'))
    sequences = {'test_0_nr':  # non-ref / target
                   'ATGATTAAGA\n',  # 9
                 }
    refSequences = {'test_0_r':  # ref / query
                    'ATGATCCAATGA\n'  # 12
                    }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0          11
    # ref        ATGATCCAATGA  query
    # exons       ****  ****
    # non ref    ATGATTAA--GA  target
    #            0          9
    #####
    pslLines = [simplePsl('+', 8, 0, 8, 10, 1, 9,
                          [4, 1, 1], [0, 4, 7], [1, 7, 8],
                          qName='ensmust0', tName='test_0_nr'),
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 1, 11, 'ensmust0', 0, '+', 1, 11,
        '128,0,0', 2, '4,4', '0,6'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 1, 9, 'ensmust0', 0, '+', 1, 9,
        '128,0,0', 2, '4,2', '0,6'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run the filter
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'indel', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))

    self.assertEqual(len(writtenTranscripts), 1)
    self.assertEqual(writtenTranscripts[0].chromosomeInterval, transcripts[0].chromosomeInterval)
    self.assertEqual(writtenTranscripts[0].exons, transcripts[0].exons)
    self.assertEqual(len(writtenTranscripts[0].annotations), 1)
    annotation = writtenTranscripts[0].annotations[0]
    self.assertEqual(len(annotation.labels), 1)
    self.assertEqual(annotation.labels[0], 'deletion')
    self.assertEqual(annotation.chromosomeInterval.start, 7)
    self.assertEqual(annotation.chromosomeInterval.stop, 9)

    self.addCleanup(removeDir, tmpDir)

  def test_indel_1(self):
    """indel should produce insertion annotations correctly."""
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('indel_1'))
    sequences =  {'test_0_r':  # non-ref / target
                  'ATGATCCAATGA\n'  # 12
                 }
    refSequences = {'test_0_nr':  # ref / query
                    'ATGATTAAGA\n',  # 9
                   }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0          9
    # ref        ATGATTAA--GA  query
    # exons       ****  ****
    # non ref    ATGATCCAATGA  target
    #            0          11
    #####
    pslLines = [simplePsl('+', 6, 0, 6, 12, 1, 11,
                          [4, 1, 1], [0, 4, 5], [1, 7, 10],
                          qName='ensmust0', tName='test_0_nr'),
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 1, 9, 'ensmust0', 0, '+', 1, 9,
        '128,0,0', 2, '4,2', '0,6'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 1, 11, 'ensmust0', 0, '+', 1, 11,
        '128,0,0', 2, '4,4', '0,6'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run the filter
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'indel', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))

    self.assertEqual(len(writtenTranscripts), 1)
    self.assertEqual(writtenTranscripts[0].chromosomeInterval, transcripts[0].chromosomeInterval)
    self.assertEqual(writtenTranscripts[0].exons, transcripts[0].exons)
    self.assertEqual(len(writtenTranscripts[0].annotations), 1)
    annotation = writtenTranscripts[0].annotations[0]
    self.assertEqual(len(annotation.labels), 1)
    self.assertEqual(annotation.labels[0], 'insertion')
    self.assertEqual(annotation.chromosomeInterval.start, 8)
    self.assertEqual(annotation.chromosomeInterval.stop, 10)

    self.addCleanup(removeDir, tmpDir)

  def test_indel_2(self):
    """indel should fix introns that are created by insertions."""
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('indel_2'))
    sequences =  {'test_0_r':  # non-ref / target
                  'ATGATCCAATGA\n'  # 12
                 }
    refSequences = {'test_0_nr':  # ref / query
                    'ATGATTAAGA\n',  # 9
                   }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0          9
    # ref        ATGATTAA--GA  query
    # exons       ****  ****
    # non ref    ATGATCCAATGA  target
    #            0          11
    #####
    pslLines = [simplePsl('+', 6, 0, 6, 12, 1, 11,
                          [4, 1, 1], [0, 4, 5], [1, 7, 10],
                          qName='ensmust0', tName='test_0_nr'),
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 1, 9, 'ensmust0', 0, '+', 1, 9,
        '128,0,0', 2, '4,2', '0,6'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 1, 11, 'ensmust0', 0, '+', 1, 11,
        '128,0,0', 3, '4,1,1', '0,6,9'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    self.assertEqual(len(transcripts), 1)
    # Add unknownSplice tags to both so that we can check that the one
    # corresponding to the removed intron gets removed but the correct
    # one stays.
    transcripts[0].annotations.append(lib_filter.TranscriptAnnotation(
      lib_filter.ChromosomeInterval(
        'test_0_nr',
        5,
        7,
        True),
      'ensmust0', ['unknownUtrSplice']))
    transcripts[0].annotations.append(lib_filter.TranscriptAnnotation(
      lib_filter.ChromosomeInterval(
        'test_0_nr',
        8,
        10,
        True),
      'ensmust0', ['unknownCdsSplice']))
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run the filter
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'indel', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))

    self.assertEqual(len(writtenTranscripts), 1)
    self.assertEqual(writtenTranscripts[0].chromosomeInterval, transcripts[0].chromosomeInterval)
    # decrease in # of exons
    self.assertEqual(len(writtenTranscripts[0].exons), 2)
    self.assertEqual([(e.start, e.stop) for e in writtenTranscripts[0].exons], [(1, 5), (7, 11)])
    # 1 insertion annotation, 1 remaining unkown splice annotation
    self.assertEqual(len(writtenTranscripts[0].annotations), 2)

    # Check the insertion annotation
    insertionAnnotations = [i for i in writtenTranscripts[0].annotations if 'insertion' in i.labels]
    self.assertEqual(len(insertionAnnotations), 1)
    insertionAnnotation = insertionAnnotations[0]
    self.assertEqual(len(insertionAnnotation.labels), 1)
    self.assertEqual(insertionAnnotation.labels[0], 'insertion')
    self.assertEqual(insertionAnnotation.chromosomeInterval.start, 8)
    self.assertEqual(insertionAnnotation.chromosomeInterval.stop, 10)

    # Check the remaining unknown splice annotation
    unknownSpliceAnnotations = [i for i in writtenTranscripts[0].annotations \
                                if 'unknownCdsSplice' in i.labels or 'unknownUtrSplice' in i.labels]
    self.assertEqual(len(unknownSpliceAnnotations), 1)
    unknownSpliceAnnotation = unknownSpliceAnnotations[0]
    self.assertEqual(len(unknownSpliceAnnotation.labels), 1)
    self.assertEqual(unknownSpliceAnnotation.labels[0], 'unknownUtrSplice')
    self.assertEqual(unknownSpliceAnnotation.chromosomeInterval.start, 5)
    self.assertEqual(unknownSpliceAnnotation.chromosomeInterval.stop, 7)

    self.addCleanup(removeDir, tmpDir)

  def test_indel_3(self):
    """indel should produce deletion annotations correctly even on exon
    boundaries."""
    makeTempDirParent()
    tmpDir = os.path.abspath(makeTempDir('indel_0'))
    sequences = {'test_0_nr':  # non-ref / target
                   'ATGATTAAGA\n',  # 9
                 }
    refSequences = {'test_0_r':  # ref / query
                    'ATGATCCAATGA\n'  # 12
                    }
    seqFile = createSequenceFile(sequences, tmpDir)
    seqDict = lib_filter.getSequences(seqFile)
    refSeqFile = createSequenceFile(refSequences, tmpDir, filename='refSeq.fa')
    refSeqDict = lib_filter.getSequences(refSeqFile)
    ##########
    #            0          11
    # ref        ATGATCCAATGA  query
    # exons       ****  ****
    # non ref    ATGATTA--AGA  target
    #            0          9
    #####
    pslLines = [simplePsl('+', 8, 0, 8, 10, 1, 9,
                          [4, 2], [0, 6], [1, 7],
                          qName='ensmust0', tName='test_0_nr'),
                ]
    pslFile = createAlignmentFile(pslLines, tmpDir)
    refTranscriptBedLines = []
    refTranscriptBedLines.append(bedLine(
        'test_0_r', 1, 11, 'ensmust0', 0, '+', 1, 11,
        '128,0,0', 2, '4,4', '0,6'))
    createBedFile(refTranscriptBedLines, 'ref.bed', tmpDir)
    transcriptBedLines = []
    transcriptBedLines.append(bedLine(
        'test_0_nr', 1, 9, 'ensmust0', 0, '+', 1, 9,
        '128,0,0', 2, '4,2', '0,6'))
    transcriptDetailsBedLines = []
    transcripts = [
      transcript for transcript in lib_filter.transcriptIterator(
        transcriptBedLines, transcriptDetailsBedLines)]
    # write transcripts to files
    outBed = os.path.join(tmpDir, 'test.bed')
    outDetailsBed = os.path.join(tmpDir, 'test_details.bed')
    testFile = lib_filter.writeTranscriptBedFile(
      transcripts, outBed)
    testDetailsFile = lib_filter.writeDetailsBedFile(
      transcripts, outDetailsBed)
    # run the filter
    with open(os.path.join(tmpDir, 'dummy.txt'), 'w') as f:
      f.write('dummy\n')
    with open(os.path.join(tmpDir, 'empty.txt'), 'w') as f:
      f.write('')
    originalGeneCheckBed = os.path.join(tmpDir, 'ref.bed')
    originalGeneCheckBedDetails = os.path.join(tmpDir, 'empty.txt')
    alignment = os.path.join(tmpDir, 'aln.psl')
    sequence = os.path.join(tmpDir, 'seq.fa')
    referenceSequence = os.path.join(tmpDir, 'refSeq.fa')
    chromSizes = os.path.join(tmpDir, 'dummy.txt')
    metaFilter.makeCall(
      'indel', 'C57B6J', 'C57B6NJ', outBed, outDetailsBed,
      originalGeneCheckBed, originalGeneCheckBedDetails,
      alignment, sequence, referenceSequence, chromSizes, tmpDir)
    # read transcripts from file
    writtenTranscripts = lib_filter.getTranscripts(
      os.path.join(tmpDir, 'out.bed'), os.path.join(tmpDir, 'out_details.bed'))

    self.assertEqual(len(writtenTranscripts), 1)
    self.assertEqual(writtenTranscripts[0].chromosomeInterval, transcripts[0].chromosomeInterval)
    self.assertEqual(writtenTranscripts[0].exons, transcripts[0].exons)
    self.assertEqual(len(writtenTranscripts[0].annotations), 1)
    annotation = writtenTranscripts[0].annotations[0]
    self.assertEqual(len(annotation.labels), 1)
    self.assertEqual(annotation.labels[0], 'deletion')
    self.assertEqual(annotation.chromosomeInterval.start, 6)
    self.assertEqual(annotation.chromosomeInterval.stop, 8)

    self.addCleanup(removeDir, tmpDir)

if __name__ == '__main__':
  unittest.main()
