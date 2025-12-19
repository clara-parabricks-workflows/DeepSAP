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

import os
import csv
import numpy as np
import pickle
import logging
import json
import logging
import argparse
from datetime import datetime

from dataclasses import dataclass, field

import torch
from torch import distributed as dist
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

import transformers
import sklearn

from typing import Optional, Dict, Sequence, Tuple, List

os.environ["TORCH_DISTRIBUTED_DEBUG"] = "OFF"
logging.disable(logging.WARNING)

os.environ["HUGGINGFACE_HOME"] = "/root/.cache/huggingface"

class SupervisedDataset(Dataset):
    """Dataset for supervised fine-tuning."""
    
    def __init__(self, 
                 data_path: str, 
                 tokenizer: transformers.PreTrainedTokenizer, 
                 kmer: int = -1):

        super(SupervisedDataset, self).__init__()

        # load data from the disk
        with open(data_path, "r") as f:
            data = list(csv.reader(f))[1:]
        if len(data[0]) == 2:
            # data is in the format of [text, label]
            logging.warning("Perform single sequence classification...")
            texts = [d[0] for d in data]
            labels = [int(d[1]) for d in data]
        elif len(data[0]) == 3:
            # data is in the format of [text1, text2, label]
            logging.warning("Perform sequence-pair classification...")
            texts = [[d[0], d[1]] for d in data]
            labels = [int(d[2]) for d in data]
        else:
            raise ValueError("Data format not supported.")
        
        # --- START OF THE ROBUST FIX ---
        if kmer != -1:
            is_distributed = dist.is_available() and dist.is_initialized()
            rank = dist.get_rank() if is_distributed else 0

            # Step 1: Only the main process generates the k-mer file.
            if rank == 0:
                logging.warning(f"Rank 0: Generating/checking {kmer}-mer file...")
                # This function should save the k-mers to a predictable file path
                load_or_generate_kmer(data_path, texts, kmer)

            # Step 2: All processes wait here until Rank 0 is done.
            # In a single-GPU run, this does nothing.
            if is_distributed:
                dist.barrier()

            # Step 3: Now that the file is guaranteed to exist, ALL processes load it.
            # This ensures `texts` is identical for every process.
            logging.warning(f"Rank {rank}: Loading {kmer}-mer as input...")
            # This function should now just load the file created by Rank 0.
            texts = load_or_generate_kmer(data_path, texts, kmer) 
        # --- END OF THE ROBUST FIX ---

        output = tokenizer(
            texts,
            return_tensors="pt",
            padding="longest",
            max_length=tokenizer.model_max_length,
            truncation=True,
        )

        self.input_ids = output["input_ids"]
        self.attention_mask = output["attention_mask"]
        self.labels = labels
        self.num_labels = len(set(labels))

    def __len__(self):
        return len(self.input_ids)

    def __getitem__(self, i) -> Dict[str, torch.Tensor]:
        return dict(input_ids=self.input_ids[i], labels=self.labels[i])

@dataclass
class DataCollatorForSupervisedDataset(object):
    """Collate examples for supervised fine-tuning."""

    tokenizer: transformers.PreTrainedTokenizer

    def __call__(self, instances: Sequence[Dict]) -> Dict[str, torch.Tensor]:
        input_ids, labels = tuple([instance[key] for instance in instances] for key in ("input_ids", "labels"))
        input_ids = torch.nn.utils.rnn.pad_sequence(
            input_ids, batch_first=True, padding_value=self.tokenizer.pad_token_id
        )
        labels = torch.Tensor(labels).long()
        return dict(
            input_ids=input_ids,
            labels=labels,
            attention_mask=input_ids.ne(self.tokenizer.pad_token_id),
        )

"""
Load or generate k-mer string for each DNA sequence. The generated k-mer string will be saved to the same directory as the original data with the same name but with a suffix of "_{k}mer".
"""
def load_or_generate_kmer(data_path: str, texts: List[str], k: int) -> List[str]:
    """Load or generate k-mer string for each DNA sequence."""
    kmer_path = data_path.replace(".csv", f"_{k}mer.json")
    if os.path.exists(kmer_path):
        logging.warning(f"Loading k-mer from {kmer_path}...")
        with open(kmer_path, "r") as f:
            kmer = json.load(f)
    else:        
        logging.warning(f"Generating k-mer...")
        kmer = [generate_kmer_str(text, k) for text in texts]
        with open(kmer_path, "w") as f:
            logging.warning(f"Saving k-mer to {kmer_path}...")
            json.dump(kmer, f)
        
    return kmer

"""
Transform a dna sequence to k-mer string
"""
def generate_kmer_str(sequence: str, k: int) -> str:
    """Generate k-mer string from DNA sequence."""
    return " ".join([sequence[i:i+k] for i in range(len(sequence) - k + 1)])

@dataclass
class ModelArguments:
    model_name_or_path: Optional[str] = field(default="facebook/opt-125m")
    use_lora: bool = field(default=False, metadata={"help": "whether to use LoRA"})
    lora_r: int = field(default=8, metadata={"help": "hidden dimension for LoRA"})
    lora_alpha: int = field(default=32, metadata={"help": "alpha for LoRA"})
    lora_dropout: float = field(default=0.05, metadata={"help": "dropout rate for LoRA"})
    lora_target_modules: str = field(default="query,value", metadata={"help": "where to perform LoRA"})


@dataclass
class DataArguments:
    data_path: str = field(default=None, metadata={"help": "Path to the training data."})
    kmer: int = field(default=-1, metadata={"help": "k-mer for input sequence. -1 means not using k-mer."})


# COPY tuned_models/transformers_cache/* /root/.cache/huggingface/transformers/

@dataclass
class TrainingArguments(transformers.TrainingArguments):
    cache_dir: Optional[str] = field(default="/root/.cache/huggingface/transformers/models--zhihan1996--DNA_bert_6/snapshots/55e0c0eb7b734c8b9b77bc083bf89eb6fbda1341/")
    run_name: str = field(default="run")
    optim: str = field(default="adamw_torch")
    model_max_length: int = field(default=512, metadata={"help": "Maximum sequence length."})
    gradient_accumulation_steps: int = field(default=1)
    per_device_train_batch_size: int = field(default=1)
    per_device_eval_batch_size: int = field(default=1)
    num_train_epochs: int = field(default=1)
    fp16: bool = field(default=False)
    logging_steps: int = field(default=100)
    save_steps: int = field(default=100)
    eval_steps: int = field(default=100)
    evaluation_strategy: str = field(default="steps")
    warmup_steps: int = field(default=50)
    weight_decay: float = field(default=0.01)
    learning_rate: float = field(default=1e-4)
    save_total_limit: int = field(default=3)
    load_best_model_at_end: bool = field(default=False)
    output_dir: str = field(default="output")
    find_unused_parameters: bool = field(default=False)
    checkpointing: bool = field(default=False)
    dataloader_pin_memory: bool = field(default=False)
    eval_and_save_results: bool = field(default=True)
    save_model: bool = field(default=False)
    seed: int = field(default=42)

"""
Manually calculate the accuracy, f1, matthews_correlation, precision, recall with sklearn.
"""
def calculate_metric_with_sklearn(logits: np.ndarray, labels: np.ndarray):
    if logits.ndim == 3:
        # Reshape logits to 2D if needed
        logits = logits.reshape(-1, logits.shape[-1])
    predictions = np.argmax(logits, axis=-1)
    valid_mask = labels != -100  # Exclude padding tokens (assuming -100 is the padding token ID)
    valid_predictions = predictions[valid_mask]
    valid_labels = labels[valid_mask]
    return {
        "accuracy": sklearn.metrics.accuracy_score(valid_labels, valid_predictions),
        "f1": sklearn.metrics.f1_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
        "matthews_correlation": sklearn.metrics.matthews_corrcoef(
            valid_labels, valid_predictions
        ),
        "precision": sklearn.metrics.precision_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
        "recall": sklearn.metrics.recall_score(
            valid_labels, valid_predictions, average="macro", zero_division=0
        ),
    }
    
"""
Compute metrics used for huggingface trainer.
""" 
def compute_metrics(eval_pred):
    logits, labels = eval_pred
    if isinstance(logits, tuple):  # Unpack logits if it's a tuple
        logits = logits[0]
    return calculate_metric_with_sklearn(logits, labels)

def predict_with_attentions(pre_trained_model_name, tuned_model_path, data_path, data_kmer, MODEL_MAX_LEN, output_dir, BATCH_SIZE=64, SEED=42, is_fp16=True, is_CUDA=True):

    os.environ["TOKENIZERS_PARALLELISM"] = "false"
    torch.cuda.empty_cache()

    # 1. Load Tokenizer (No change)
    tokenizer = transformers.AutoTokenizer.from_pretrained(
        pre_trained_model_name,
        model_max_length=MODEL_MAX_LEN,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True,
    )

    # 2. Load Dataset and Collator (No change)
    predict_dataset = SupervisedDataset(tokenizer=tokenizer,
                                        data_path=os.path.join(data_path, "dev.csv"),
                                        kmer=data_kmer)
    data_collator = DataCollatorForSupervisedDataset(tokenizer=tokenizer)

    # 3. Load Model
    #    - Set output_attentions=True here for clarity.
    print("Loading model...")
    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        tuned_model_path,
        num_labels=predict_dataset.num_labels,
        trust_remote_code=True,
        output_attentions=True, # Set to True for consistency
        output_hidden_states=True,
    )

    if is_CUDA:
        model.to('cuda')
    model.eval() # Set the model to evaluation mode

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # 4. Manual Prediction Loop (Cleaned up)
    all_attention_scores = []
    all_hidden_states    = []
    all_tokens           = []
    all_probs            = [] # <-- NEW: List to store prediction probabilities

    dataloader = DataLoader(predict_dataset, batch_size=BATCH_SIZE, collate_fn=data_collator)

    print("Running predictions...")
    for i, batch in enumerate(dataloader):
        # NOTE: This only processes the first sequence in each batch.
        # If batch size > 1, the tokens will not match the full attention/hidden states.
        tokens = tokenizer.convert_ids_to_tokens(batch['input_ids'][0].tolist())
        all_tokens.append(tokens)
        
        # Move the whole batch to the GPU
        inputs = {k: v.to('cuda') for k, v in batch.items()} if is_CUDA else batch

        # Forward pass to get logits, hidden states, and attention scores
        with torch.no_grad():
            outputs = model(**inputs)

        # --- EXTRACT AND STORE PROBABILITIES ---
        # Get logits and apply softmax to convert to probabilities
        probs = torch.softmax(outputs.logits, dim=-1)
        # Move to CPU and append to list
        all_probs.append(probs.cpu())
        # ----------------------------------------

        # Detach and move tensors to CPU before appending to prevent VRAM leak
        all_attention_scores.append(tuple(att.cpu() for att in outputs.attentions))
        all_hidden_states.append(tuple(hs.cpu() for hs in outputs.hidden_states))

        # This break is useful for debugging, keeping it.
        if i > 1000:
            print("Stopping after 1000 batches for debugging.")
            break
            
    # 5. Save all outputs

    # --- SAVE PREDICTION PROBABILITIES ---
    # Concatenate list of tensors into one big tensor
    final_probs = torch.cat(all_probs, dim=0).numpy()
    output_probs_file = os.path.join(output_dir, "predictions_probs.npy")
    print(f"**** Writing prediction probabilities into {output_probs_file}")
    np.save(output_probs_file, final_probs)
    # -----------------------------------

    output_states_file = os.path.join(output_dir, "hidden_states.pkl")
    print(f"**** Writing hidden states into {output_states_file}")
    with open(output_states_file, 'wb') as handle:
        pickle.dump(all_hidden_states, handle, protocol=pickle.HIGHEST_PROTOCOL)

    output_attn_file = os.path.join(output_dir, "attns.pkl")
    print(f"**** Writing attentions into {output_attn_file}")
    with open(output_attn_file, 'wb') as handle:
        pickle.dump(all_attention_scores, handle, protocol=pickle.HIGHEST_PROTOCOL)

    output_tokens_file = os.path.join(output_dir, "tokens.pkl")
    print(f"**** Writing tokens into {output_tokens_file}")
    with open(output_tokens_file, 'wb') as handle:
        pickle.dump(all_tokens, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print("Prediction and artifact extraction complete.")

def predict(pre_trained_model_name, tuned_model_path, data_path , data_kmer, MODEL_MAX_LEN, output_dir, BATCH_SIZE=64, SEED=42, is_fp16=True):
    """
    Runs inference using a fine-tuned transformer model.
    """
    # os.environ["TOKENIZERS_PARALLELISM"] = "false"
    torch.cuda.empty_cache()

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tLoading tokenizer...")
    tokenizer = transformers.AutoTokenizer.from_pretrained(
        "/DNA_bert_6/", 
        model_max_length=MODEL_MAX_LEN,
        padding_side="right",
        use_fast=True,
        trust_remote_code=False,
        local_files_only=True
    )

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tLoading inference dataset...")
    predict_dataset = SupervisedDataset(tokenizer=tokenizer,
                                        data_path=os.path.join(data_path, "dev.csv"),
                                        kmer=data_kmer)
   
    data_collator = DataCollatorForSupervisedDataset(tokenizer=tokenizer)
 
    # print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tLoading model from '{tuned_model_path}'...")
    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        tuned_model_path,
        num_labels=predict_dataset.num_labels,
        output_attentions=False,
        trust_remote_code=False,
        local_files_only=True
    )

    prediction_args = transformers.TrainingArguments(
        output_dir=output_dir,                      # Directory where all outputs (predictions, checkpoints) will be saved.
        per_device_eval_batch_size=BATCH_SIZE,      # Batch size for evaluation on each available GPU or CPU.
        save_strategy="steps",                      # Checkpoint saving strategy; not critical for a prediction-only script.
        evaluation_strategy="no",                   # Disables automatic evaluation during a training run.
        dataloader_drop_last=False,                 # Ensures the last, smaller batch is processed instead of being dropped.
        overwrite_output_dir=True,                  # Allows the script to write over files in an existing output directory.
        fp16=is_fp16,                               # Enables 16-bit floating point (half-precision) for faster inference.
        remove_unused_columns=True,                 # Automatically removes data columns not used by the model's forward pass.
        seed=SEED,                                  # Sets the random seed for reproducibility.
        report_to="none",                           # Suppress wandb/tensorboard warnings
    )
    
    trainer = transformers.Trainer(model=model,
                                   tokenizer=tokenizer,
                                   args=prediction_args,
                                   compute_metrics=compute_metrics,
                                   train_dataset=None,
                                   eval_dataset=None,
                                   data_collator=data_collator)

    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tStarting inference ...")
    prediction_output = trainer.predict(predict_dataset, metric_key_prefix="predict_")
    
    logits = prediction_output.predictions

    # --- DEBUGGING BLOCK ---
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[DEBUG]\tShape of logits from trainer: {logits.shape}")
    if not isinstance(logits, np.ndarray) or logits.ndim != 2:
        raise ValueError(f"Expected a 2D numpy array of logits from the trainer, but received an object of type {type(logits)} with shape {getattr(logits, 'shape', 'N/A')}.")
    # --- END DEBUGGING BLOCK ---

    softmax = torch.nn.Softmax(dim=1)
    probs = softmax(torch.tensor(logits, dtype=torch.float32)).numpy()
    
    output_predict_file = os.path.join(prediction_args.output_dir, "probs")
    np.save(output_predict_file, probs)
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tWrote prediction probabilities to {output_predict_file}.npy")

    # print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tStarting inference ...")
    # predictions, _, _ = trainer.predict(predict_dataset, metric_key_prefix="predict_")

    # softmax = torch.nn.Softmax(dim=1)
    
    # try:
    #     probs_array = np.array(predictions[0])
    # except:
    #     probs_array = np.array(predictions)
        
    # probs = softmax(torch.tensor(probs_array, dtype=torch.float32)).numpy()
    
    # output_predict_file = os.path.join(prediction_args.output_dir, "probs")
    # np.save(output_predict_file, probs)
    # print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}]\t[LOG]\tWrote prediction probabilities to {output_predict_file}.npy")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Tuning Nucleotide Transformers', usage="python tune_NT.py -m nucleotide-transformer-2.5b-multi-species -d /SpliceUp/datasets/RefSeq400 -o /SpliceUp/tune_models/NT_2.5B_850g_multi_species__Splice__RefSeq400 -x 100 -n 10 -r 6.1e-5 -b 8")
    parser.add_argument('-m','--model',  help='Pre-trained model name',                 required=True)
    parser.add_argument('-t','--tuned',  help='Tuned model path',                       required=True)
    parser.add_argument('-d','--data',   help='Path to data folder',                    required=True)
    parser.add_argument('-k','--kmer',   help='Kmer size',                              required=True, default=-1)
    parser.add_argument('-x','--max',    help='Max sequence length after tokenization', required=True)
    parser.add_argument('-o','--out',    help='Path to output folder',                  required=True)
    # parser.add_argument('-l','--LoRA',   help='Using LoRA',                             required=False, default=False, type=bool)
    parser.add_argument('-s','--seed',   help='Randomness seed',                        required=False, default=42)
    
    parser.add_argument('--fp16',           help='Use fp16', action='store_true')
    parser.add_argument('--no-fp16',        help="Don't use fp16", dest='fp16', action='store_false')
    parser.add_argument('--process',        help='Initiate process', action='store_true')
    parser.add_argument('--no-process',     help="Don't initiate process", dest='process', action='store_false')

    parser.set_defaults(fp16=True)
    parser.set_defaults(process=True)

    args = parser.parse_args()
    
    """
    usage: 
    
    python ../scripts/predict.py \
        --model "InstaDeepAI/nucleotide-transformer-2.5b-multi-species" \
        --tuned "/Users/fberakdar/Library/Mobile Documents/com~apple~CloudDocs/Nvidia/code/DeepSNAP/tuned_models/NT_2.5B_850g_multi_species__Splice__RefSeq400/" \
        --data "/Users/fberakdar/Library/Mobile Documents/com~apple~CloudDocs/Nvidia/code/DeepSNAP/datasets/RefSeq400_dummy/" \
        --kmer -1 \
        --max 100 \
        --batch 32 \
        --out ${output} &> ${output}/log_predict.txt
    """

    predict(args.model, args.tuned, args.data , int(args.kmer), int(args.max), args.out, int(args.batch), int(args.seed), bool(args.fp16))