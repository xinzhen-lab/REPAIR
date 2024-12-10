import subprocess

def run_training():
    dataset = "brats_hgg_lgg"
    miss_rate = 0.4
    command = ["python", "demo.py", "--dataset", dataset, "--miss_rate", str(miss_rate)]

    try:
        # 执行命令
        subprocess.run(command, check=True)
        print("Training completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running the command: {e}")

if __name__ == "__main__":
    run_training()
