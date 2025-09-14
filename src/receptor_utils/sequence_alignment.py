# Functions for creating formatted sequence alignments

# Copyright (c) 2021 William Lees

# This source code, and any executable file compiled or derived from it, is governed by the European Union Public License v. 1.2,
# the English version of which is available here: https://perma.cc/DK5U-NDVE

from receptor_utils import simple_bio_seq as simple


def create_alignment(sequences, sequence_type='V', codon_wrap=20, v_coords=None, output_format='text'):
    """Create a formatted sequence alignment showing variations from the first sequence.
    
    :param sequences: dict of sequences in the format returned by read_fasta
    :type sequences: dict
    :param sequence_type: type of sequence (V, D, or J) - currently only V is implemented
    :type sequence_type: str
    :param codon_wrap: number of codons per line before wrapping
    :type codon_wrap: int
    :param v_coords: coordinates for V-gene regions as (cdr1_start, cdr1_end, cdr2_start, cdr2_end, cdr3_start)
    :type v_coords: tuple or None
    :param output_format: output format ('text' or 'html' - currently only text is implemented)
    :type output_format: str
    :return: formatted alignment string
    :rtype: str
    """
    
    if sequence_type not in ['V', 'D', 'J']:
        raise ValueError("sequence_type must be 'V', 'D', or 'J'")
    
    if output_format != 'text':
        raise NotImplementedError("Only text output format is currently implemented")
    
    # Set default V coordinates if not provided (only used for V sequences)
    if v_coords is None and sequence_type == 'V':
        v_coords = (27, 38, 56, 65, 105)  # Default IMGT V-gene coordinates (1-based)
    # Validate coordinates only for V sequences
    if sequence_type == 'V':
        if len(v_coords) != 5:
            raise ValueError("v_coords must contain exactly 5 integers")
        
        # Validate coordinates are ascending
        for i in range(1, len(v_coords)):
            if v_coords[i] <= v_coords[i-1]:
                raise ValueError("v_coords values must be in ascending order")
    
    if not sequences:
        return ""
    
    # Sort sequences alphabetically by name
    sorted_seq_names = sorted(sequences.keys())
    
    # Get the first sequence as reference
    reference_name = sorted_seq_names[0]
    reference_seq = sequences[reference_name]
    
    # Calculate dynamic starting position based on longest sequence name
    max_name_length = max(len(name) for name in sorted_seq_names)
    nucleotide_start_pos = max_name_length + 5  # 5 spaces after longest name
    
    # Normalize gaps - convert both . and - to .
    normalized_sequences = {}
    for name, seq in sequences.items():
        normalized_sequences[name] = seq
    
    reference_seq = normalized_sequences[reference_name]
    
    # Translate reference sequence to amino acids
    reference_aa = simple.translate(reference_seq)
    
    # Create alignment output
    alignment_lines = []
    
    # Calculate total codons needed
    total_codons = len(reference_seq) // 3
    
    # Process alignment in chunks based on codon_wrap
    for chunk_start in range(0, total_codons, codon_wrap):
        chunk_end = min(chunk_start + codon_wrap, total_codons)
        
        # Add header lines for this chunk
        if sequence_type == 'V':
            alignment_lines.extend(_create_chunk_header(chunk_start + 1, chunk_end, v_coords, nucleotide_start_pos))
        else:
            alignment_lines.extend(_create_chunk_header_no_cdr(chunk_start + 1, chunk_end, nucleotide_start_pos))
        
        # Add reference sequence line
        ref_line = _create_reference_line(reference_name, reference_seq, reference_aa,
                                          chunk_start, chunk_end, nucleotide_start_pos, max_name_length, sequence_type)
        alignment_lines.append(ref_line)
        alignment_lines.append("")  # Empty line after reference
        
        # Add comparison sequences
        for seq_name in sorted_seq_names[1:]:
            if seq_name in normalized_sequences:
                comp_line = _create_comparison_line(seq_name, normalized_sequences[seq_name],
                                                  reference_seq, reference_aa,
                                                  chunk_start, chunk_end, nucleotide_start_pos, max_name_length, sequence_type)
                alignment_lines.append(comp_line)
                alignment_lines.append("")  # Empty line after each sequence
        
        alignment_lines.append("")  # Extra empty line between chunks
    
    return '\n'.join(alignment_lines)


def _create_chunk_header(start_codon, end_codon, v_coords, nucleotide_start_pos):
    """Create header lines for a chunk of the alignment."""
    lines = []
    
    # Add CDR region labels - positioned above codon numbers
    region_line = " " * (nucleotide_start_pos + 4*(end_codon - start_codon + 1))  # Initial spacing
    
    # Handle CDR1
    if v_coords[0] <= end_codon and v_coords[1] >= start_codon:  # CDR1 overlaps this chunk
        cdr1_start_in_chunk = max(v_coords[0], start_codon)
        cdr1_end_in_chunk = min(v_coords[1], end_codon)
        cdr1_start_pos = nucleotide_start_pos + (cdr1_start_in_chunk - start_codon) * 4
        cdr1_span_in_chunk = cdr1_end_in_chunk - cdr1_start_in_chunk + 1
        cdr1_label_length = (cdr1_span_in_chunk - 1) * 4 + 3
        
        # Create appropriate label for this chunk portion
        if v_coords[0] >= start_codon and v_coords[1] <= end_codon:
            # Complete CDR1 in this chunk
            underscores_each_side = (cdr1_label_length - len("CDR1")) // 2
            label = "_" * underscores_each_side + "CDR1" + "_" * underscores_each_side
            if len(label) < cdr1_label_length:
                label += "_"
        elif v_coords[0] >= start_codon:
            # CDR1 starts in this chunk but continues beyond
            label = "_" * (cdr1_label_length - len("CDR1")) + "CDR1"
        elif v_coords[1] <= end_codon:
            # CDR1 started before this chunk but ends here
            label = "CDR1" + "_" * (cdr1_label_length - len("CDR1"))
        else:
            # CDR1 spans entire chunk (started before, continues after)
            label = "_" * cdr1_label_length
        
        region_line = region_line[:cdr1_start_pos] + label + region_line[cdr1_start_pos + len(label):]
    
    # Handle CDR2
    if v_coords[2] <= end_codon and v_coords[3] >= start_codon:  # CDR2 overlaps this chunk
        cdr2_start_in_chunk = max(v_coords[2], start_codon)
        cdr2_end_in_chunk = min(v_coords[3], end_codon)
        cdr2_start_pos = nucleotide_start_pos + (cdr2_start_in_chunk - start_codon) * 4
        cdr2_span_in_chunk = cdr2_end_in_chunk - cdr2_start_in_chunk + 1
        cdr2_label_length = (cdr2_span_in_chunk - 1) * 4 + 3
        
        # Create appropriate label for this chunk portion
        if v_coords[2] >= start_codon and v_coords[3] <= end_codon:
            # Complete CDR2 in this chunk
            underscores_each_side = (cdr2_label_length - len("CDR2")) // 2
            label = "_" * underscores_each_side + "CDR2" + "_" * underscores_each_side
            if len(label) < cdr2_label_length:
                label += "_"
        elif v_coords[2] >= start_codon:
            # CDR2 starts in this chunk but continues beyond
            label = "_" * (cdr2_label_length - len("CDR2")) + "CDR2"
        elif v_coords[3] <= end_codon:
            # CDR2 started before this chunk but ends here 
            label = "CDR2" + "_" * (cdr2_label_length - len("CDR2"))
        else:
            # CDR2 spans entire chunk (started before, continues after)
            label = "_" * cdr2_label_length
        
        region_line = region_line[:cdr2_start_pos] + label + region_line[cdr2_start_pos + len(label):]
    
    # Handle CDR3
    if v_coords[4] <= end_codon:  # CDR3 starts in or before this chunk
        cdr3_start_in_chunk = max(v_coords[4], start_codon)
        cdr3_start_pos = nucleotide_start_pos + (cdr3_start_in_chunk - start_codon) * 4
        
        # CDR3 extends to end of sequence, so just mark the start
        if v_coords[4] >= start_codon:
            label = "__CDR3_"
            region_line = region_line[:cdr3_start_pos] + label + region_line[cdr3_start_pos + len(label):]
    
    lines.append(region_line)
    lines.append("")  # Empty line
    
    # Add codon numbers - center above middle nucleotide of each codon
    # Sequence names are formatted dynamically, nucleotides start at nucleotide_start_pos
    # Middle of first codon is at nucleotide_start_pos + 1, subsequent codons are 4 chars apart
    number_line = " " * (nucleotide_start_pos + 1)  # Position to align with middle nucleotide of first codon
    for i in range(start_codon, end_codon + 1):
        # Calculate spacing to center number above middle nucleotide
        num_str = str(i)
        # Each codon is 4 chars wide, center the number
        if len(num_str) == 1:
            spacing = "   "  # 3 spaces after single digit
        elif len(num_str) == 2:
            spacing = "  "   # 2 spaces after double digit
        else:
            number_line = number_line[:-1]
            spacing = "  "    # 1 spaces after triple digit
        number_line += num_str + spacing
    
    lines.append(number_line)
    
    return lines


def _create_chunk_header_no_cdr(start_codon, end_codon, nucleotide_start_pos):
    """Create header lines for a chunk of D or J alignment (no CDR regions)."""
    lines = []
    
    # Add empty line where CDR regions would be (for consistent spacing)
    lines.append("")
    lines.append("")  # Empty line
    
    # Add codon numbers - center above middle nucleotide of each codon
    # Sequence names are formatted dynamically, nucleotides start at nucleotide_start_pos
    # Middle of first codon is at nucleotide_start_pos + 1, subsequent codons are 4 chars apart
    number_line = " " * (nucleotide_start_pos + 1)  # Position to align with middle nucleotide of first codon
    for i in range(start_codon, end_codon + 1):
        # Calculate spacing to center number above middle nucleotide
        num_str = str(i)
        # Each codon is 4 chars wide, center the number
        if len(num_str) == 1:
            spacing = "   "  # 3 spaces after single digit
        elif len(num_str) == 2:
            spacing = "  "   # 2 spaces after double digit
        else:
            spacing = " "    # 1 space after triple digit
        number_line += num_str + spacing
    
    lines.append(number_line)
    
    return lines


def _create_reference_line(seq_name, nucleotide_seq, amino_acid_seq, start_codon, end_codon, nucleotide_start_pos, max_name_length, sequence_type='V'):
    """Create the reference sequence line showing both nucleotides and amino acids."""
    # Extract the relevant portion of the sequence
    start_nt = start_codon * 3
    end_nt = end_codon * 3
    
    # Handle gaps in the nucleotide sequence
    nt_ungapped = nucleotide_seq
    aa_ungapped = amino_acid_seq
    
    # Create amino acid line first - only for V and J sequences, not D sequences
    if sequence_type != 'D':
        # Position to center above middle nucleotide of each codon
        # Sequence names are formatted dynamically, nucleotides start at nucleotide_start_pos
        # Middle of first codon is at nucleotide_start_pos + 1, subsequent amino acids are 4 chars apart
        aa_line = " " * (nucleotide_start_pos + 1)  # Position to align with middle nucleotide of first codon
        chunk_aa = aa_ungapped[start_codon:end_codon]
        for i, aa in enumerate(chunk_aa):
            aa_line += f"{aa}   "
    
    # Create the nucleotide line with sequence name
    line = f"{seq_name:<{max_name_length}}     "  # Name padded to max length + 5 spaces
    
    # Add nucleotides for this chunk
    chunk_nt = nt_ungapped[start_nt:end_nt]
    
    # Format nucleotides in codons
    nt_part = ""
    for i in range(0, len(chunk_nt), 3):
        codon = chunk_nt[i:i+3]
        if len(codon) == 3:
            nt_part += codon.upper() + " "
        elif len(codon) > 0:
            nt_part += codon.upper().ljust(3) + " "
    
    # Return with or without amino acid line depending on sequence type
    if sequence_type == 'D':
        return line + nt_part
    else:
        return aa_line + "\n" + line + nt_part


def _create_comparison_line(seq_name, nucleotide_seq, ref_nt_seq, ref_aa_seq, start_codon, end_codon, nucleotide_start_pos, max_name_length, sequence_type='V'):
    """Create a comparison line showing only differences from reference."""
    # Extract relevant portions
    start_nt = start_codon * 3
    end_nt = end_codon * 3
    
    ref_nt_ungapped = ref_nt_seq
    comp_nt_ungapped = nucleotide_seq
    
    # Get chunks
    ref_chunk_nt = ref_nt_ungapped[start_nt:end_nt]
    comp_chunk_nt = comp_nt_ungapped[start_nt:end_nt] if len(comp_nt_ungapped) > start_nt else ""
    
    # Create amino acid differences line - only for V and J sequences, not D sequences
    if sequence_type != 'D':
        aa_line = " " * (nucleotide_start_pos + 1)  # Position to align with middle nucleotide of first codon
        aa_changes = []
    
    # Create the nucleotide line
    line = f"{seq_name:<{max_name_length}}     "  # Name padded to max length + 5 spaces
    
    # Compare codons
    nt_part = ""
    
    for i in range(0, len(ref_chunk_nt), 3):
        ref_codon = ref_chunk_nt[i:i+3] if i+3 <= len(ref_chunk_nt) else ref_chunk_nt[i:]
        comp_codon = comp_chunk_nt[i:i+3] if i+3 <= len(comp_chunk_nt) else comp_chunk_nt[i:] if i < len(comp_chunk_nt) else ""
        
        if len(ref_codon) == 3 and len(comp_codon) == 3:
            if ref_codon == comp_codon:
                nt_part += "--- " if (comp_codon != '...' and comp_codon != '---') else '... '
                if sequence_type != 'D':
                    aa_changes.append("")
            else:
                if sequence_type == 'D':
                    # For D sequences, just show nucleotide changes without amino acid analysis
                    codon_display = ""
                    for j in range(3):
                        if ref_codon[j] == comp_codon[j]:
                            codon_display += "-"
                        else:
                            codon_display += comp_codon[j].lower()
                    nt_part += codon_display + " "
                else:
                    # For V and J sequences, check if it's a silent mutation
                    ref_aa = simple.translate(ref_codon)
                    comp_aa = simple.translate(comp_codon)
                    
                    if ref_aa == comp_aa:
                        # Silent mutation - show only nucleotide changes
                        codon_display = ""
                        for j in range(3):
                            if ref_codon[j] == comp_codon[j]:
                                codon_display += "-"
                            else:
                                codon_display += comp_codon[j].lower()
                        nt_part += codon_display + " "
                        aa_changes.append("")
                    else:
                        # Non-silent mutation - show nucleotide changes and AA
                        codon_display = ""
                        for j in range(3):
                            if ref_codon[j] == comp_codon[j]:
                                codon_display += "-"
                            else:
                                codon_display += comp_codon[j].lower()
                        nt_part += codon_display + " "
                        aa_changes.append(comp_aa)
        elif len(comp_codon) == 0:
            # Sequence ends here
            break
        else:
            # Partial codon
            nt_part += comp_codon.ljust(3) + " " if comp_codon else ""
            if sequence_type != 'D':
                aa_changes.append("")
    
    # Handle output based on sequence type
    if sequence_type == 'D':
        # For D sequences, only return nucleotide line
        return line + nt_part
    else:
        # For V and J sequences, add amino acid changes to the aa_line
        for i, aa_change in enumerate(aa_changes):
            if aa_change:
                aa_line += f"{aa_change}   "
            else:
                aa_line += "    "
        
        # Create the full output
        if any(aa_changes):
            return aa_line + "\n" + line + nt_part
        else:
            return " " * (nucleotide_start_pos + 1) + "\n" + line + nt_part
