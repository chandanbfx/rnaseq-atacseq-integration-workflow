# Common functions and variables
import os

def get_log(sample, step, type="log"):
    return os.path.join(config["log_dir"], f"{sample}_{step}.{type}")

def ensure_dir(d):
    os.makedirs(d, exist_ok=True)
