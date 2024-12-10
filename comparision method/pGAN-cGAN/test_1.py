import subprocess
import torch

print(torch.__version__)

def train():
    # 设置训练相关的参数
    dataroot = "datasets/brats/"
    name = "cGAN_run"
    model = "cGAN"
    lambda_A = 100
    lambda_B = 100
    batch_size = 1
    output_nc = 1
    input_nc = 1  # 根据具体任务调整
    gpu_ids = 0
    niter = 50
    niter_decay = 50
    save_epoch_freq = 25
    checkpoints_dir = "checkpoints/brats_t2_flair/"
    dataset_mode = "unaligned_mat"

    # 训练命令
    train_command = [
        "python", "cGAN.py",
        "--dataroot", dataroot,
        "--name", name,
        "--model", model,
        "--output_nc", str(output_nc),
        "--input_nc", str(input_nc),
        "--gpu_ids", str(gpu_ids),
        "--niter", str(niter),
        "--niter_decay", str(niter_decay),
        "--save_epoch_freq", str(save_epoch_freq),
        "--lambda_A", str(lambda_A),
        "--lambda_B", str(lambda_B),
        "--checkpoints_dir", checkpoints_dir,
        "--dataset_mode", dataset_mode,
        "--training"
    ]

    print("开始训练...")
    subprocess.run(train_command)

def test():
    # 设置测试相关的参数
    dataroot = "datasets/REMBRANDT/"
    name = "cGAN_run"
    output_nc = 1
    input_nc = 1
    how_many = 50000
    results_dir = "result/REMBRANDT/"
    checkpoints_dir = "checkpoints/brats_t1_t1ce/"

    # 测试命令
    test_command = [
        "python", "cGAN.py",
        "--dataroot", dataroot,
        "--name", name,
        "--model", "cGAN",
        "--phase", "test",
        "--output_nc", str(output_nc),
        "--input_nc", str(input_nc),
        "--how_many", str(how_many),
        "--results_dir", results_dir,
        "--checkpoints_dir", checkpoints_dir
    ]

    print("开始测试...")
    subprocess.run(test_command)

if __name__ == "__main__":
    # train()  # 调用训练函数
    test()   # 调用测试函数
    print("训练和测试完成。")
