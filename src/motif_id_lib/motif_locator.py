import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

degenerate_map = {
    'A': ['A'], 
    'T': ['T'], 
    'G': ['G'], 
    'C': ['C'],
    'N': ['A','T','G','C'],
    'R': ['A','G'],
    'Y': ['C','T'],
    'W': ['A','T'],
    'S': ['G','C'],
    'K': ['G','T'],
    'M': ['A','C'],
    'B': ['C','G','T'],
    'D': ['A','G','T'], 
    'H': ['A','C','T'], 
    'V': ['A','C','G']
}

class Motifs:
  def __init__(self, enzymes):
    self.enzymes_dict = enzymes
    self.plasmid = None
    self.motif_results = {}

  def array_set(self, plasmid: str):
    '''
    Purpose
    -------
    Transforms the plasmid sequence into a Numpy array.

    Args
    ----
    plasmid sequence as a string

    '''
    self.plasmid = np.array(list(plasmid))

  def motif_search(self, motif: str) -> np.ndarray:
    '''
    Purpose
    -------
    Finds the motif locations in the plasmid sequence.

    Args
    ----
    motif sequence as a string

    Returns
    -------
    Numpy array containing the starting indices of all motif matches in the plasmid sequence.
    
    Notes and Example:
    - Step 1: Translate the motif sequence into a list containing sets for each base. 
    - Step 2: Create sliding windows over the plasmid sequence. The window size is the length of the motif.
    - Step 3: Initialize match mask for each window in the motif (all True).
    - Step 4: Filter the mask windows for each base in the motif. The match mask is updated to false for any window that does not match the motif at that postiion.
    - Step 5: Return the starting indices of all "True" matches.

    For the motif "AR" and the plasmid sequence "AGTGCAT":
    - Step 1: The motif "AR" is translated into the list [{A}, {A, G}] based on the degenerate map.
    - Step 2: Sliding windows of size 2 are created over the plasmid sequence, resulting in windows ["AG", "GT", "TG", "GC", "CA", "AT"].
    - Step 3: A match mask for each window is initialized as [True, True, True, True, True, True]. 
    - Step 4: The match mask is updated based on the motif:
      Checks for the first base 'A' in the motif, the mask is updated to [True, False, False, False, True, True] because only "AG" and "AT" have 'A' at the first position.
      Then, checks for the second base 'R' (which can be 'A' or 'G'), the mask is further updated to [True, False, False, False, False, False] because "AG" has 'G' at the second position, and "AT" has 'T' which does not match 'R'.
    - Step 5: The starting indices of "True" matches are returned, which in this case would be index 0 corresponding to the match "AG" in the plasmid sequence.
    '''

    # Step 1
    updated_motif = [set(degenerate_map[base]) for base in motif]
    
    # Step 2
    windows = sliding_window_view(self.plasmid, len(motif))
    
    # Step 3
    matches = np.ones(len(windows), dtype=bool)
    
    # Step 4
    for i, motif_bases in enumerate(updated_motif):
        matches &= np.isin(windows[:, i], list(motif_bases))
    
    # Step 5
    positions = np.where(matches)[0]
    return positions

  def get_motif_results(self):
    '''
    Purpose
    -------
    Create a dictionary containing the enzyme names as keys and the values as a list of the motif, the cut site, the count of motif observations, and the locations of the motif in the plasmid sequence.

    Returns
    -------
    Dictionary with enzyme names as keys and the values as a list of the motif, cut site, observation counts, and locations in the plasmid sequence. 
    '''
    for enzyme in self.enzymes_dict:
        motif = self.enzymes_dict[enzyme][0]
        cut_site = self.enzymes_dict[enzyme][1]
        positions = self.motif_search(motif)
        count = len(positions)
        self.motif_results[enzyme] = [motif, cut_site, count, positions.tolist()]
    
    return self.motif_results
