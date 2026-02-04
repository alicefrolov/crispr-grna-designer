### **grna_designer.py:**
```python
#!/usr/bin/env python3
"""
CRISPR gRNA Designer
Tool for identifying and evaluating potential CRISPR guide RNAs.
Author: Alice Frolov
"""

def find_pam_sites(sequence):
    """Find all PAM (NGG) sites in the sequence."""
    pam_sites = []
    sequence = sequence.upper()
    
    for i in range(len(sequence) - 2):
        # Look for NGG pattern (N = any nucleotide)
        if sequence[i+1:i+3] == 'GG':
            pam_sites.append(i)
    
    return pam_sites

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def check_homopolymers(sequence, max_length=4):
    """Check for homopolymer runs that may cause off-target effects."""
    for base in ['A', 'T', 'G', 'C']:
        if base * max_length in sequence:
            return True
    return False

def evaluate_grna(grna_seq):
    """Evaluate gRNA quality based on multiple criteria."""
    gc_content = calculate_gc_content(grna_seq)
    has_homopolymer = check_homopolymers(grna_seq)
    
    score = 0
    warnings = []
    
    # GC content scoring (optimal: 40-60%)
    if 40 <= gc_content <= 60:
        score += 3
    elif 30 <= gc_content <= 70:
        score += 2
        warnings.append("GC content suboptimal")
    else:
        score += 1
        warnings.append("GC content outside recommended range")
    
    # Homopolymer check
    if has_homopolymer:
        score -= 1
        warnings.append("Contains homopolymer run (potential off-target risk)")
    else:
        score += 2
    
    # Check for specific problematic sequences
    if 'TTTT' in grna_seq:
        warnings.append("Contains TTTT (Pol III terminator)")
        score -= 2
    
    return score, warnings

def design_grnas(target_sequence, grna_length=20):
    """Design gRNA candidates from a target sequence."""
    target_sequence = target_sequence.upper()
    pam_sites = find_pam_sites(target_sequence)
    
    candidates = []
    
    for pam_pos in pam_sites:
        # Extract 20bp upstream of PAM
        start = pam_pos - grna_length
        if start >= 0:
            grna_seq = target_sequence[start:pam_pos]
            
            if len(grna_seq) == grna_length:
                pam = target_sequence[pam_pos:pam_pos+3]
                gc_content = calculate_gc_content(grna_seq)
                score, warnings = evaluate_grna(grna_seq)
                
                candidates.append({
                    'sequence': grna_seq,
                    'pam': pam,
                    'position': start,
                    'gc_content': gc_content,
                    'score': score,
                    'warnings': warnings
                })
    
    # Sort by score (best first)
    candidates.sort(key=lambda x: x['score'], reverse=True)
    
    return candidates

def print_results(candidates):
    """Print gRNA candidates in a formatted table."""
    if not candidates:
        print("No suitable gRNA candidates found.")
        return
    
    print(f"\nFound {len(candidates)} gRNA candidate(s)\n")
    print("="*80)
    
    for i, candidate in enumerate(candidates[:10], 1):  # Show top 10
        print(f"\nCandidate #{i}")
        print("-"*80)
        print(f"gRNA Sequence:  5'-{candidate['sequence']}-3'")
        print(f"PAM:            {candidate['pam']}")
        print(f"Position:       {candidate['position']}")
        print(f"GC Content:     {candidate['gc_content']:.1f}%")
        print(f"Quality Score:  {candidate['score']}/5")
        
        if candidate['warnings']:
            print(f"Warnings:       {', '.join(candidate['warnings'])}")
        else:
            print(f"Warnings:       None (Good candidate!)")
    
    print("\n" + "="*80)
    print(f"\nTop candidate: 5'-{candidates[0]['sequence']}-3'")
    print(f"Recommended for experimental validation.")

def main():
    """Main function to run the gRNA designer."""
    import sys
    
    print("CRISPR gRNA Designer")
    print("="*80)
    
    if len(sys.argv) < 2:
        print("\nUsage: python grna_designer.py <target_sequence>")
        print('Example: python grna_designer.py "ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG"')
        sys.exit(1)
    
    target_sequence = sys.argv[1]
    
    print(f"\nTarget Sequence: {target_sequence}")
    print(f"Length: {len(target_sequence)} bp\n")
    
    candidates = design_grnas(target_sequence)
    print_results(candidates)

if __name__ == "__main__":
    main()