import os
import sys
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import yaml
import pandas as pd
import re
from tqdm import tqdm
from abnumber import Chain
from abnumber.exceptions import ChainParseError
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqIO import parse
from Bio.Align import substitution_matrices, PairwiseAligner
from pymol import cmd


def get_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    return aligner


def get_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    return aligner


def get_chain_type(seq, antigen_ref, antigen_min_score, antigen_min_score_per_residue):
    chain_info = {}
    resi_map = {}

    try:
        # Create a Chain object with the cleaned sequence
        ab_chain = Chain(seq, assign_germline=True, scheme="imgt")

        # Align the sequence to the IMGT scheme
        std_seq = ["-"] * 129  # Assuming IMGT alignment length is 129 positions
        for i, (position, aa) in enumerate(ab_chain.positions.items()):
            # Extract the numeric part of the IMGT position
            key = re.sub('[a-zA-Z]', '', str(position))
            if key.isdigit():
                imgt_pos = int(key)
                std_seq[imgt_pos - 1] = aa
                resi_map[i + 1] = imgt_pos

        imgt_aligned_seq = "".join(std_seq)

        chain_info.update(
            {
                "type": f"{'heavy' if ab_chain.chain_type == 'H' else 'light'}",
                "species": ab_chain.species,
                "v_gene": ab_chain.v_gene,
                "j_gene": ab_chain.j_gene,
                "seq": imgt_aligned_seq,
                "resi_map": resi_map,
            }
        )
    except ChainParseError:
        aligner = get_aligner()
        aln = aligner.align(antigen_ref, seq)

        if len(aln) == 0:
            return None

        aln = aln[0]

        ag_start = aln.aligned[0][0][0]
        ag_end = aln.aligned[0][-1][1]

        seq_aligned_len = aln.aligned[1][-1][1] - aln.aligned[1][0][0] + 1

        if (
            aln.score < antigen_min_score
            or aln.score / seq_aligned_len < antigen_min_score_per_residue
        ):
            return None

        resi_map = {}
        for ref_seg, seq_seg in zip(aln.aligned[0], aln.aligned[1]):
            assert ref_seg[1] - ref_seg[0] == seq_seg[1] - seq_seg[0]

            for j in range(seq_seg[0], seq_seg[1]):
                resi_map[str(j + 1)] = str(ref_seg[0] + 1 + j - seq_seg[0])

        chain_info.update(
            {
                "type": "antigen",
                "seq": seq[aln.aligned[1][0][0] : (aln.aligned[1][-1][1] + 1)],
                "start": ag_start,
                "end": ag_end,
                "score": aln.score,
                "score_per_residue": format(aln.score / seq_aligned_len, ".2f"),
                "resi_map": resi_map,
            }
        )
    except AssertionError:
        return None

    return chain_info


def get_chains(
    pdb_dict: dict,
    ref: str,
    antigen_min_score: int,
    antigen_min_score_per_residue: float,
):
    resi_nums = {}
    for cid, strand, auth_num, num, ins_code in zip(
        pdb_dict["_pdbx_poly_seq_scheme.entity_id"],
        pdb_dict["_pdbx_poly_seq_scheme.pdb_strand_id"],
        pdb_dict["_pdbx_poly_seq_scheme.pdb_seq_num"],
        pdb_dict["_pdbx_poly_seq_scheme.seq_id"],
        pdb_dict["_pdbx_poly_seq_scheme.pdb_ins_code"],
    ):
        # Only extract the numeric part of the residue number
        if strand not in resi_nums:
            resi_nums[strand] = {}
        auth_num_numeric = ''.join(filter(str.isdigit, auth_num))
        if auth_num_numeric:  # Add only if numeric part exists
            resi_nums[strand][auth_num_numeric] = num



    poly_chains = pdb_dict["_entity_poly.entity_id"]
    poly_desc = {}
    for i in poly_chains:
        for j, cid in enumerate(pdb_dict["_entity.id"]):
            if cid == i and cid not in poly_desc:
                poly_desc[cid] = pdb_dict["_entity.pdbx_description"][j]
    poly_strand_id = pdb_dict["_entity_poly.pdbx_strand_id"]
    poly_seq = pdb_dict["_entity_poly.pdbx_seq_one_letter_code_can"]
    poly_type = pdb_dict["_entity_poly.type"]

    ab_chains = {}
    ag_chains = {}

    for cid, strand, seq, chain_type in zip(
        poly_chains, poly_strand_id, poly_seq, poly_type
    ):
        if chain_type != "polypeptide(L)":
            continue
        for chars in ["\n"]:
            seq = seq.replace(chars, "")

        if len(seq) == 0:
            continue

        chain_res = get_chain_type(
            seq,
            antigen_ref=ref,
            antigen_min_score=antigen_min_score,
            antigen_min_score_per_residue=antigen_min_score_per_residue,
        )
        if chain_res is None:
            continue
        chain_res["description"] = poly_desc[cid]

        for strand_single in strand.split(","):
            chain_info = chain_res.copy()
            if chain_info["type"] == "heavy" or chain_info["type"] == "light":
                chain_info["resi_num"] = resi_nums[strand_single]  # Store residue numbering for antibodies
                ab_chains[strand_single] = chain_info

            elif chain_info["type"] == "antigen":
                chain_info["resi_num"] = resi_nums[strand_single]
                ag_chains[strand_single] = chain_info

    return ab_chains, ag_chains


def get_paratope(ab_chains, ag_chains):
    """
    Calculate the paratope sites on the antibody chain.
    """
    results = {}
    cmd.remove('solvent')  # Remove solvent molecules

    for ab_chain in ab_chains:
        cmd.select(
            "paratope", f"chain {ab_chain} within 3.5 of chain {','.join(ag_chains)}"
        )

        _sites = set([x.resi for x in cmd.get_model("paratope").atom])
        
        # Keep only numeric parts of the residue identifiers
        _sites = [resi for resi in _sites if resi.isdigit()]

        # Convert the residue identifiers to integers for further use
        _sites = list(map(int, _sites))

        if len(_sites) > 0:
            results[ab_chain] = _sites

    cmd.delete("all")
    return results


def get_epitope(
    entry, ab_chain, ab_chains, ag_chains, thres_ag, thres_ab,ab_range=None,resi_num=None, imgt_numbering=None):
    cmd.load(entry)
    cmd.remove('solvent')
    results = {}
    contact_ab_chains = []

    # find contacting glycans and assign to corresponding antigen chains

    cmd.select("antibody", f"c. {ab_chain}" + ("" if ab_range is None else f" and resi {ab_range[0]}-{ab_range[1]}"))
    cmd.select("epitope_glycan", f"organic within {thres_ag} of antibody")
    glycan_resi = set([(x.chain, x.resi) for x in cmd.get_model("epitope_glycan").atom])
    glycan_asn_resi = {}

    for chain, resi in glycan_resi:
        _n_ext = 0
        cmd.select("_ext_glycan", f"resi {resi} and c. {chain}")
        while len(cmd.get_model("name ND2 and _ext_glycan").atom) == 0 and _n_ext < 10:
            cmd.select("_ext_glycan", "_ext_glycan xt. 2")
            _n_ext += 1

        if _n_ext < 10:
            _atoms = cmd.get_model("name ND2 and _ext_glycan").atom
            for _atom in _atoms:
                chain_asn, resi_asn = _atom.chain, _atom.resi
                if chain_asn in ag_chains:
                    if chain_asn not in glycan_asn_resi:
                        glycan_asn_resi[chain_asn] = set()
                    glycan_asn_resi[chain_asn].add(resi_asn)

    # find epitope sites on each antigen chain
    for ag_chain in ag_chains:
        cmd.select(
            "epitope",
            f"polymer.protein and c. {ag_chain} within {thres_ag} of antibody",
        )

        _sites = set([x.resi for x in cmd.get_model("epitope").atom])
        _sites |= glycan_asn_resi.get(ag_chain, set())
        _sites = list(_sites)

        if len(_sites) > 0:
            results[ag_chain] = _sites

    # find pairing antibody chain
    for ab_c in ab_chains:
        cmd.select(
            "epitope",
            f"polymer.protein and c. {ab_chain} within {thres_ab} of c. {ab_c}",
        )
        if len(cmd.get_model("epitope").atom) > 0 and ab_c != ab_chain:
            contact_ab_chains.append((ab_c, len(cmd.get_model("epitope").atom)))

    if len(contact_ab_chains) > 0:
        contact_ab_chains = sorted(contact_ab_chains, key=lambda x: x[1], reverse=True)[0][0]
    else:
        contact_ab_chains = ""

    paratopes = get_paratope(ab_chain, ag_chains)


    cmd.delete("all")
    return results, contact_ab_chains,paratopes


def get_antigen_type(antigens, epitope_sites):
    counts = {}
    for ag, (start, end) in antigens.items():
        for site in epitope_sites:
            if site[0] != "<" and start <= int(site) <= end:
                if ag in counts:
                    counts[ag] += 1
                else:
                    counts[ag] = 1

    return "+".join(counts.keys()) if len(counts) > 0 else "unknown"


def parse_cif(
    entry: Path,
    reference: str,
    antigens: dict[str, tuple[int, int]],
    antigen_contact_threshold: float,
    antibody_pair_contact_threshold: float,
    antigen_min_score: int,
    antigen_min_score_per_residue: float,
):
    pdb_id = entry.stem
    pdb_dict = dict(MMCIF2Dict(entry).items())
    ab_chains, ag_chains = get_chains(
        pdb_dict,
        ref=reference,
        antigen_min_score=antigen_min_score,
        antigen_min_score_per_residue=antigen_min_score_per_residue,
    )

    if len(ab_chains) == 0 or len(ag_chains) == 0:
        return pdb_id

    results = []
    involved_antigen_strands = set()
    for ab, abinfo in ab_chains.items():
        imgt_numbering = abinfo.get("imgt_numbering", None)
        resi_num = abinfo.get("resi_num", None)
        epitopes, contact_ab_chain, paratopes= get_epitope(
            entry,
            ab,
            [_c for _c, _props in ab_chains.items() if _props.get('type') != abinfo.get('type')],
            ag_chains.keys(),
            thres_ag=antigen_contact_threshold,
            thres_ab=antibody_pair_contact_threshold,
            resi_num=resi_num,
            imgt_numbering=imgt_numbering
            # ab_range=(1, len(abinfo["seq"])),
        )

        _antigen_types = []
        _strands = []
        _epitope_sites = []
        _paratope_sites = []

        for strand in sorted(epitopes.keys()):
            # print(ag_chains[strand])
            resi_map = ag_chains[strand]["resi_map"]
            resi_num = ag_chains[strand]["resi_num"]

            numbers = sorted(
                [(resi_num[x], x) for x in epitopes[strand]], key=lambda x: int(x[0])
            )

            epitope_sites = [
                resi_map[num_x] if num_x in resi_map else f"<{x}>"
                for num_x, x in numbers
            ]
            _antigen_types.append(get_antigen_type(antigens, epitope_sites))
            _epitope_sites.append("+".join(epitope_sites))
            _strands.append(strand)

            involved_antigen_strands.add(strand)
        
        _paratope_sites = []
        
        _paratope_sites = []

        _paratope_sites = []

        for strand in paratopes:
            print(f"Processing strand: {strand}")  # Debug print for the current strand
            
            # Get residue numbering information
            resi_map = ab_chains[strand].get("resi_map", {})
            resi_num = ab_chains[strand].get("resi_num", {})
            resi_map = {int(k): v for k, v in resi_map.items()}

           
            # Check the contents of paratopes and resi_num for debugging
            print(f"paratopes for strand {strand}: {paratopes.get(strand, [])}")
            print(f"resi_num for strand {strand}: {resi_map}")

            # Normalize paratope site IDs
            paratope_sites_raw = paratopes.get(strand, [])
            
            # Convert `paratope_sites_raw` to strings to match with `resi_num` keys
            paratope_sites_raw_int = [int(site) for site in paratope_sites_raw]            
            print(f"paratope_sites_raw_str for strand {strand}: {paratope_sites_raw_int}")
            # Only consider valid residue numbers that exist in `resi_map`
            matched_keys = []
            for key in paratope_sites_raw_int:
                print(f"key: {key}")
                if key in resi_map:
                    matched_keys.append(resi_map[key])
                    print(f"matched key: {resi_map[key]}")
                
            if matched_keys!= []:
                sorted_matched_keys = sorted(matched_keys, key=lambda x: int(x))

                paratope_sites = [str(pos) for pos in sorted_matched_keys]

                # Join sites with "+" separator to match the required format
                _paratope_sites.append("+".join(paratope_sites))
            print(f"_paratope_sites: {_paratope_sites}")

            
            
            
        report_features = [
            "start",
            "end",
            "seq",
            "score",
            "score_per_residue",
            "description",
        ]
        _ag_info = {
            "antigen_strand": "|".join(_strands),
            "epitope_sites": "|".join(_epitope_sites),
            "paratope_sites_true": "|".join(_paratope_sites),
            "antigen_type": "|".join(_antigen_types),
            **{
                f"antigen_{y}": "|".join([str(ag_chains[x][y]) for x in _strands])
                for y in report_features
            },
        }

        results.append(
            {
                "pdb_id": pdb_id,
                "strand": ab,
                "contact_ab_chain": contact_ab_chain,
                **_ag_info,
                **abinfo,
            }
        )
    
    if len(results) == 0 or len(involved_antigen_strands) == 0:
        return pdb_id
    return pd.DataFrame(results)


parser = argparse.ArgumentParser(description="Process PDB mmCIF files")
parser.add_argument("--config", "-c", type=str, default="config.yaml", help="Configuration yaml file")
parser.add_argument("--input_dir", "-i", type=str, default=None, help="Input directory containing PDB mmCIF files")
parser.add_argument("--output_dir", "-o", type=str, default=None, help="Output directory")
parser.add_argument("--reference", "-r", type=str, default=None, help="Reference sequence fasta file")
parser.add_argument("--ncpu", "-n", type=int, default=None, help="Number of CPUs to use")


def main(args):
    with open(args.config, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
    # set working directory to the parent of config file

    for k, v in vars(args).items():
        if v is not None and k not in config and k != "config":
            config[k] = v
    
    HIV_ENV_REF = next(parse(config["reference"], "fasta")).seq
    print(
        f"Reference HIV Env Sequence:\n{str(HIV_ENV_REF)}\nLength: {len(HIV_ENV_REF)}",
        file=sys.stdout,
        flush=True,
    )

    df = []
    failed_pdb_entries = []

    pdbs = list(Path(config["input_dir"]).glob("*.cif"))

    print(f"Processing {len(pdbs)} PDB entries", file=sys.stdout, flush=True)

    num_proc = min(int(config["ncpu"]), len(pdbs), os.cpu_count())
    if num_proc <= 0:
        num_proc = os.cpu_count()

    if num_proc > 1:
        with ProcessPoolExecutor(max_workers=config["ncpu"]) as executor:
            futures = [
                executor.submit(
                    parse_cif,
                    **{
                        "entry": cif,
                        "reference": HIV_ENV_REF,
                        "antigens": config["antigens"],
                        "antigen_contact_threshold": config[
                            "antigen_contact_threshold"
                        ],
                        "antibody_pair_contact_threshold": config[
                            "antibody_pair_contact_threshold"
                        ],
                        "antigen_min_score": config["antigen_min_score"],
                        "antigen_min_score_per_residue": config[
                            "antigen_min_score_per_residue"
                        ],
                    },
                )
                for cif in pdbs
            ]
            for future in tqdm(futures):
                x = future.result(timeout=20)
                if isinstance(x, str):
                    failed_pdb_entries.append(x)
                else:
                    df.append(x)
    else:
        for cif in tqdm(pdbs):
            print(f"Processing {cif}", file=sys.stdout, flush=True)
            x = parse_cif(
                cif,
                HIV_ENV_REF,
                config["antigens"],
                config["antigen_contact_threshold"],
                config["antibody_pair_contact_threshold"],
                config["antigen_min_score"],
                config["antigen_min_score_per_residue"],
            )
            if isinstance(x, str):
                failed_pdb_entries.append(x)
            else:
                df.append(x)

    df = pd.concat(df).sort_values(["pdb_id", "description"])


    output_dir = Path(config["output_dir"])
    output_dir.mkdir(exist_ok=True, parents=True)
    df.to_csv(output_dir / "processed.csv", index=None)

    with open(output_dir / "failed_pdb_entries.txt", "w", encoding="utf-8") as f:
        for entry in failed_pdb_entries:
            f.write(entry + "\n")


if __name__ == "__main__":
    arguments = parser.parse_args()
    main(arguments)
