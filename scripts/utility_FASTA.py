# Copyright 2025 NVIDIA CORPORATION
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import numpy as np
import subprocess
import gzip
from datetime import datetime
from pyfaidx import Fasta
import torch
import transformers
import os, sys, re
import contextlib

def parse_fasta(path:str):
    """Parsing FASTA file

    Args:
        path (string): Path to FASTA file

    Returns:
        sequences (dictionary): dictionary of seauence name to string of bases
    """
    
    print("[{}]\t[LOG]\tParsing FASTA file '{}'".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), path))
    
    sequences = dict()
    buffer = ""
    header = ""

    file   = open(path, "r")
    zipped = False
    
    if path[-3:] == ".gz":
        file = gzip.open(path, 'rb')
        zipped = True
        
    while True:
        line = file.readline()
        
        if len(line) == 0:
            break
        
        if zipped:
            line = line.decode()
        
        if line[0] == '>':
            if len(buffer) != 0:
                sequences[header] = buffer
            
            buffer = ""
            if ' ' in line:
                header = line.split(' ')[0][1:]
            else:
                header = line[1:-1]
        else:
            buffer += line[:-1]
    
    if len(buffer) != 0:
        sequences[header] = buffer
                
    file.close()
    return sequences 




# def get_reverse_complement(dna_sequence):
#     """Returns the reverse complement of a DNA sequence."""
#     complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
#     return "".join(complement[base] for base in reversed(dna_sequence))

# def encode_probabilities_to_byte(probs):
#     """
#     Encodes a 4-element probability vector into a single byte by binning
#     each of the 4 probabilities into a 2-bit value.
#     """
#     binned_values = []
#     for p in probs:
#         if p < 0.25:
#             binned_values.append(0b00)
#         elif p < 0.50:
#             binned_values.append(0b01)
#         elif p < 0.75:
#             binned_values.append(0b10)
#         else:
#             binned_values.append(0b11)

#     final_byte_int = (binned_values[0] << 6) | (binned_values[1] << 4) | (binned_values[2] << 2) | binned_values[3]
#     return final_byte_int.to_bytes(1, 'big')

# def process_and_write_batch(sequences_to_predict, predictions_batch, tsv_writer, bin_writer, trainer, predict_args):
#     """
#     Predicts a batch of sequences using the pre-loaded trainer and writes results.
#     This version operates entirely in memory, avoiding temporary files.
#     """

#     class SupervisedDataset(torch.utils.data.Dataset):
#         """
#         The Dataset class accepts a list of sequences directly in memory,
#         avoiding the need to read from a temporary file.
#         """
#         def __init__(self, tokenizer, sequences):
#             self.sequences = sequences
#             self.tokenizer = tokenizer
#             self.num_labels = 3

#         def __len__(self):
#             return len(self.sequences)

#         def __getitem__(self, idx):
#             return self.tokenizer(self.sequences[idx], return_tensors="pt", padding="max_length", truncation=True)
        
#     if not sequences_to_predict:
#         print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tNo valid sequences to process in this chunk.")
#     else:
#         print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tPredicting on {len(sequences_to_predict)} sequences in memory...")
        
#         # Prepare sequences for the model directly in memory
#         all_sequences = [item[2] for item in sequences_to_predict]
#         kmer_sequences = [kmer_split(seq, predict_args['kmer']) for seq in all_sequences]
        
#         # Add dummy sequences for models trained on 3 labels
#         kmer_sequences.extend(["N N N N N", "N N N N N", "N N N N N"])

#         # Create a dataset directly from the in-memory list
#         batch_dataset = SupervisedDataset(tokenizer=trainer.tokenizer, sequences=kmer_sequences)
        
#         # Run prediction
#         prediction_output = trainer.predict(batch_dataset, metric_key_prefix="predict_")
#         logits = prediction_output.predictions

#         softmax = torch.nn.Softmax(dim=1)
#         probs = softmax(torch.tensor(logits, dtype=torch.float32)).numpy()
        
#         # Update the main batch array with results, excluding dummy predictions
#         probs = probs[:-3] 
#         for i, (batch_idx, strand_type, _) in enumerate(sequences_to_predict):
#             if strand_type == 'forward':
#                 predictions_batch[batch_idx, 0:2] = probs[i, 0:2]
#             else: # reverse
#                 predictions_batch[batch_idx, 2:4] = probs[i, 0:2]

#     # Write the completed batch to files
#     encoded_bytes = []
#     for i in range(len(predictions_batch)):
#         preds = predictions_batch[i]
#         tsv_writer.write(f"{preds[0]:.6f}\t{preds[1]:.6f}\t{preds[2]:.6f}\t{preds[3]:.6f}\n")
#         encoded_bytes.append(encode_probabilities_to_byte(preds))
    
#     bin_writer.write(b"".join(encoded_bytes))


# ========= Fast DNA helpers =========
def kmer_split(sequence: str, kmer: int) -> str:
    """Return space-separated overlapping k-mers for DNABERT-style tokenization."""
    if kmer <= 0 or len(sequence) < kmer:
        return ""
    return " ".join(sequence[i:i+kmer] for i in range(len(sequence) - kmer + 1))

_RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def fast_revcomp(s: str) -> str:
    return s.translate(_RC_TABLE)[::-1]

def _make_base_lut():
    lut = np.full(256, -1, dtype=np.int8)
    lut[ord('A')] = 0; lut[ord('C')] = 1; lut[ord('G')] = 2; lut[ord('T')] = 3
    return lut

_BASE_LUT = _make_base_lut()

def _dinuc_hits(seq_upper: str, motif: str) -> np.ndarray:
    """Return 0-based positions i where seq[i:i+2] == motif (uppercase)."""
    b = np.frombuffer(seq_upper.encode('ascii'), dtype=np.uint8)
    if b.size < 2:
        return np.empty(0, dtype=np.int64)
    a0 = _BASE_LUT[b[:-1]]; a1 = _BASE_LUT[b[1:]]
    valid = (a0 >= 0) & (a1 >= 0)
    code = (a0 << 2) + a1
    tgt = (_BASE_LUT[ord(motif[0])] << 2) + _BASE_LUT[ord(motif[1])]
    hits = valid & (code == tgt)
    return np.flatnonzero(hits)

def _rc_2mer(d: str) -> str:
    return fast_revcomp(d)

def parse_motifs_arg(s: str) -> list[str]:
    """Accept 'GTAG ATAC GCAG' or 'GTAG,GCAG,ATAC' -> ['GTAG','ATAC','GCAG'] (deduped, validated)."""
    if not s:
        return []
    parts = [p.strip().upper() for p in re.split(r"[,\s]+", s) if p.strip()]
    seen, out = set(), []
    for p in parts:
        if len(p) == 4 and all(c in "ACGT" for c in p) and p not in seen:
            seen.add(p); out.append(p)
    return out

def parse_motif_pairs(motif_pairs: list[str]):
    """
    e.g. ["GTAG","ATAC","GCAG"] -> donors={'GT','AT','GC'}, acceptors={'AG','AC'}
    (first 2 for donor, last 2 for acceptor)
    """
    donors = set()
    acceptors = set()
    for m in motif_pairs:
        m = m.strip().upper()
        if len(m) != 4 or any(c not in "ACGT" for c in m):
            raise ValueError(f"Bad motif '{m}'. Use 4-mers like GTAG/ATAC/GCAG.")
        donors.add(m[:2]); acceptors.add(m[2:])
    return donors, acceptors

# ========= Precision / speed helpers =========
def _map_dtype(name: str | None):
    name = (name or 'fp32').lower()
    if name == 'fp16': return torch.float16
    if name == 'bf16': return torch.bfloat16
    return torch.float32

@contextlib.contextmanager
def _maybe_autocast(enabled: bool, dtype: torch.dtype | None = None):
    if enabled:
        with torch.cuda.amp.autocast(dtype=dtype if dtype in (torch.float16, torch.bfloat16) else None):
            yield
    else:
        yield

# ========= Model inference =========
def _predict_class_probs(texts, tokenizer, model, batch_size, class_idx, predict_args):
    """
    Returns np.array [len(texts)] of probabilities for class_idx.
    Fixed-length windows => default to no padding.
    """
    if not texts:
        return np.empty(0, dtype=np.float32)

    # tokenizer options (for completeness; with fixed windows use padding='none')
    padding_mode = (predict_args.get('padding', 'none') or 'none').lower()
    if padding_mode == 'none':
        padding_kw = False
        truncation_kw = False
        max_length = None
    elif padding_mode == 'longest':
        padding_kw = True
        truncation_kw = False
        max_length = None
    else:  # 'max_length'
        padding_kw = 'max_length'
        truncation_kw = True
        max_length = int(predict_args.get('max_seq', 0)) or None

    add_special = bool(predict_args.get('add_special_tokens', True))

    # precision/amp
    amp_enabled = bool(predict_args.get('amp', True))
    amp_dtype = _map_dtype(predict_args.get('dtype', 'fp32'))

    device = next(model.parameters()).device
    probs_cls = np.empty((len(texts),), dtype=np.float32)
    model.eval()

    mb = max(1, int(batch_size))
    w = 0
    with torch.inference_mode(), _maybe_autocast(amp_enabled, amp_dtype):
        for s in range(0, len(texts), mb):
            chunk = texts[s:s+mb]

            enc = tokenizer(
                chunk,
                padding=padding_kw,
                truncation=truncation_kw,
                max_length=max_length,
                add_special_tokens=add_special,
                return_tensors="pt"
            )

            input_ids = enc["input_ids"].to(device, non_blocking=True)
            attention_mask = enc["attention_mask"].to(device, non_blocking=True)

            logits = model(input_ids=input_ids, attention_mask=attention_mask).logits
            p = torch.softmax(logits, dim=1)[:, class_idx].float().cpu().numpy()
            probs_cls[w:w+len(p)] = p
            w += len(p)
    return probs_cls

# ========= Main: scan candidates & write per-chromosome CSV =========
def predict_splice_junctions_from_FASTA(
    fasta_path,
    window_size,
    n_datapoints_per_run,   # candidate batch size
    chr_coverage,
    output_path,
    predict_args
):
    """
    Produces per-chromosome CSV.GZ with columns:
      pos,motif,class,strand,probability

    Conventions (forward reference index i is first base of matched 2-mer):
      Donor + (GT/AT/GC): pos = i
      Donor - (AC/AT/GC): pos = i+1
      Acceptor + (AG/AC): pos = i+1
      Acceptor - (CT/GT): pos = i
    """
    os.makedirs(output_path, exist_ok=True)

    # --- Global matmul speed knobs (TF32) ---
    if predict_args.get('allow_tf32', True):
        torch.backends.cuda.matmul.allow_tf32 = True
        try:
            torch.set_float32_matmul_precision("high")
        except Exception:
            pass

    # Window flanks (asymmetric allowed)
    donor_left  = int(predict_args.get("donor_left", 73))
    accept_left = int(predict_args.get("accept_left", 75))
    donor_right  = window_size - donor_left  - 2
    accept_right = window_size - accept_left - 2
    if donor_right < 0 or accept_right < 0:
        raise ValueError("window_size too small for chosen left flanks (need >= left + 2).")

    # Label indices for classifier head
    donor_label_idx    = int(predict_args.get("donor_label_idx", 0))
    acceptor_label_idx = int(predict_args.get("acceptor_label_idx", 1))

    batch_size = int(predict_args['batch'])
    kmer       = int(predict_args['kmer'])
    min_prob   = float(predict_args.get('min_prob', 0.10))
    build_rc   = "flank-swapped"

    # Motif sets
    donors_fwd     = set(predict_args['donors_fwd'])
    acceptors_fwd  = set(predict_args['acceptors_fwd'])
    donors_rc      = { _rc_2mer(d) for d in donors_fwd }      # e.g., GT->AC
    acceptors_rc   = { _rc_2mer(a) for a in acceptors_fwd }   # e.g., AG->CT

    # Load tokenizer/model once
    local_only  = bool(predict_args.get('local_files_only', True))
    load_dtype  = _map_dtype(predict_args.get('dtype', 'fp32'))
    low_cpu     = bool(predict_args.get('low_cpu_mem_usage', True))
    device_map  = predict_args.get('device_map', 'cuda')

    tokenizer = transformers.AutoTokenizer.from_pretrained(
        predict_args['model_name'],
        use_fast=True,
        trust_remote_code=predict_args.get('trust_remote_code', False),
        local_files_only=local_only,
    )
    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        predict_args['model_path'],
        num_labels=int(predict_args.get('num_labels', 3)),
        torch_dtype=load_dtype,
        low_cpu_mem_usage=low_cpu,
        trust_remote_code=predict_args.get('trust_remote_code', False),
        local_files_only=local_only,
        device_map=None if device_map in (None, 'none') else {'': 'cuda'} if device_map == 'cuda' else device_map,
    )
    if device_map in (None, 'none'):
        model.to('cuda' if torch.cuda.is_available() else 'cpu')

    # ---- Force PyTorch SDPA fast attention (flash/mem) ----
    import torch.backends.cuda as tbc
    try:
        tbc.sdp_kernel(enable_flash=True, enable_mem_efficient=True, enable_math=False)
        # print("[ATTN] SDPA flash/mem enabled")
    except Exception as e:
        print(f"[ATTN] Could not set SDPA kernels: {e}")

    # ---- BetterTransformer fastpath (uses SDPA under the hood) ----
    # if predict_args.get('bettertransformer', True):
    #     try:
    #         model = model.to_bettertransformer()
    #         print("[ATTN] BetterTransformer enabled")
    #     except Exception as e:
    #         print(f"[ATTN] BetterTransformer not available: {e}")

    # Optional: environment report
    # def report_backends(model):
    #     print(f"[VERS] torch={torch.__version__} cuda={torch.version.cuda}")
    #     try:
    #         import triton; print(f"[VERS] triton={triton.__version__}")
    #     except Exception as e:
    #         print(f"[VERS] triton not importable: {e}")
    #     print("[ATTN] impl flag =", getattr(model.config, "_attn_implementation", None))
    #     print(f"[CFG] build_rc={build_rc}")
    # report_backends(model)

    genome = Fasta(fasta_path)
    # print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[CFG]\t {}")
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[CFG]\tDonors+={sorted(donors_fwd)} Donors-={sorted(donors_rc)}  Acceptors+={sorted(acceptors_fwd)} Acceptors-={sorted(acceptors_rc)}")

    for chrom_name in genome.keys():
        seq = str(genome[chrom_name]).upper()
        chrom_len = len(seq)
        max_pos = min(chrom_len, int(chr_coverage))
        out_csv = os.path.join(output_path, f"{chrom_name}.csv.gz")

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[CHR]\t{chrom_name} len={chrom_len:,} → {out_csv}")
        with gzip.open(out_csv, "wt") as f:
            # Header (no chromosome column as files are per-chromosome)
            f.write("pos,motif,class,strand,probability\n")

            # ----- helper: scan one class (D or A) without duplicates -----
            def _collect_rows_for_class(class_tag: str, fwd_2mers: set[str], rc_2mers: set[str],
                                        left_flank: int, right_flank: int, label_idx: int):
                """
                We scan the union of forward motifs and rc-only motifs to avoid duplicates.
                Palindromes (e.g., GC, AT) live in fwd_2mers and will be treated once as 'fwd'.

                For reverse hits, RC window building depends on build_rc:
                  - 'naive'         : slice [i-left : i+2+right] then RC
                  - 'flank-swapped' : slice [i-right : i+2+left] then RC (to re-center motif at 'left', which is exactlty how the model was tuned.)
                """
                rc_only = rc_2mers - fwd_2mers
                scan_motifs = list(sorted(fwd_2mers | rc_only))  # deterministic

                hits = []
                for m in scan_motifs:
                    idxs = _dinuc_hits(seq[:max_pos], m)
                    if idxs.size:
                        kind = "fwd" if m in fwd_2mers else "rc"
                        hits.append((kind, m, idxs))

                for kind, m, idxs in hits:
                    for s in range(0, idxs.size, n_datapoints_per_run):
                        block = idxs[s:s+n_datapoints_per_run]

                        meta = []      # (pos_report, motif2, class_tag, kind)
                        texts = []     # windows to feed (already oriented as the model will see)

                        # Class-specific baseline flanks
                        L_fwd = donor_left if class_tag == 'D' else accept_left
                        R_fwd = donor_right if class_tag == 'D' else accept_right

                        for i in block.astype(int, copy=False):
                            # Reported position logic
                            if class_tag == 'D':
                                pos_report = i if kind == "fwd" else (i + 1)
                            else:  # 'A'
                                pos_report = (i + 1) if kind == "fwd" else i

                            # Choose slicing flanks depending on strand and build_rc
                            if kind == "fwd":
                                L, R = L_fwd, R_fwd
                            else:
                                if build_rc == "flank-swapped":
                                    L, R = R_fwd, L_fwd  # swap on rc
                                else:
                                    L, R = L_fwd, R_fwd  # naive

                            start = i - L
                            end   = i + 2 + R  # exclusive
                            if start < 0 or end > max_pos:
                                continue
                            w = seq[start:end]
                            if 'N' in w:
                                continue

                            if kind == "fwd":
                                w_feed = w
                                expected2 = m
                                left_idx = L
                            else:
                                w_feed = fast_revcomp(w)
                                expected2 = _rc_2mer(m)
                                left_idx = L  # after swap or not, we check at index L

                            meta.append((pos_report, m, class_tag, kind))
                            texts.append(kmer_split(w_feed, kmer))

                        if not meta:
                            continue

                        probs = _predict_class_probs(texts, tokenizer, model, batch_size, label_idx, predict_args)

                        for (pos_report, motif2, klass, knd), p in zip(meta, probs):
                            strand = '+' if knd == 'fwd' else '-'
                            if p >= min_prob:
                                f.write(f"{pos_report},{motif2},{klass},{strand},{p:.6f}\n")

            # Donors (class='D')
            _collect_rows_for_class('D', donors_fwd, donors_rc, donor_left, donor_right, donor_label_idx)
            # Acceptors (class='A')
            _collect_rows_for_class('A', acceptors_fwd, acceptors_rc, accept_left, accept_right, acceptor_label_idx)

        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[DONE]\t{chrom_name}")