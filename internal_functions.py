import sys
import random
import collections
import subprocess
import matplotlib.pyplot as plt
import msprime
import numpy as np
import dataclasses
import tskit
import os
import re
import argparse
import math
import pandas as pd
import shutil as sh
from Bio import SeqIO
import glob
import seaborn as sns
import scipy.stats as stats
from Bio.Seq import Seq
from Bio import AlignIO
import Bio.Align
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceCalculator


# Create a function to check the data type of inputs
#NOTE: NOT WORKING.

def check_dat(variable_value, variable_name, expected_type):
    """Checks the type of a variable and prints a message."""
    if isinstance(variable_value, expected_type):
        print(f"Variable '{variable_name}' is of the expected data type: {expected_type.__name__}")
    else:
        print(f"Variable '{variable_name}' is NOT the expected data type. Found: {type(variable_value).__name__}")
        sys.exit()

