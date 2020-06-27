# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''A dumb wrapper for RNAfold from the Vienna RNA Package for Python 3.
Uses subprocess.Popen to open a stdin/stdout/stderr stream to RNAfold, with
the --noPS option, allowing inline use of RNAfold on RNA strings with output
wrapped in a simple object to separate sequence, energy and structure.
'''
import collections
import subprocess
import random
# Warn if Vienna not found.
from shutil import which
if not which("RNAfold"):
    from sys import stderr
    print("WARNING: PySplicer could not locate the RNAfold binary from the ViennaRNA package. Without ViennaRNA installed, structural issues cannot be diagnosed or resolved, and DNA/RNA instability may result. To put it bluntly, your designs may suck. Install ViennaRNA for best results (or any results at all)!", file=stderr)
del(which)

# Used for testing:
randseq = lambda n: ''.join([random.choice("ACGU") for x in range(0,n)])

RNAStructure = collections.namedtuple("RNAStructure",["structure","energy"])

class RNAFoldError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAfold.'
    pass

class RNAFoldOutput:
    'Wraps the two-line output from RNAfold and extracts sequence, structure and energy.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        self.sequence = output_lines[0]
        structure = output_lines[1].split(None,1)[0].strip()
        energy = float(output_lines[1].rsplit("(",1)[1].strip("()").strip())
        self.folding = RNAStructure(structure, energy)

def RNAfold(sequence, *args):
    # Note that RNAfold auto-converts "T" to "U" so this is unnecessary in Python
    # This behaviour can be overridden if calculation of DNA is desired.
    rnaf = subprocess.Popen(["RNAfold","--noPS"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNAFoldErr(folderr)
    output = RNAFoldOutput(foldout)
    #print("Debug: Energy was", output.folding.energy,"- RNAfold output was:\n",foldout) # Debug
    return output
    
class RNASuboptError(BaseException):
    'Used to wrap and raise messages from stderr when calling RNAsubopt.'
    pass

class RNASuboptOutput:
    'Wraps the muti-line output from RNAsubopt and extracts sequence, structures and energies.'
    def __init__(self, rnafold_output):
        output_lines = rnafold_output.strip().splitlines()
        # Top line of RNAsubopt has sequence and two numbers; split, take first.
        self.sequence = output_lines.pop(0).strip().split()[0]
        self.foldings = []
        for structure in output_lines:
            # Create a namedtuple of structure/energy by splitting line.
            structure, energy = structure.strip().split()
            self.foldings.append([structure, energy])

def RNAsubopt(sequence, *args):
    rnaf = subprocess.Popen(["RNAsubopt"] + list(args),
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            # Universal Newlines effectively allows string IO.
                            universal_newlines=True)        
    foldout, folderr = rnaf.communicate(sequence)
    if folderr:
        raise RNASuboptError(folderr)
    return RNASuboptOutput(foldout)
    
    