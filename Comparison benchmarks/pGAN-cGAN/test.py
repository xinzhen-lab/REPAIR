import subprocess
import torch

print(torch.__version__)



def train():
    # 设置训练相关的参数
    dataroot = "datasets/brats/"
    name = "pGAN_run"
    which_direction = "AtoB"
    lambda_A = 100
    batch_size = 1
    output_nc = 1
    input_nc = 3#层数
    gpu_ids = 0
    niter =50
    niter_decay = 50
    save_epoch_freq = 25
    lambda_vgg = 100
    checkpoints_dir = "checkpoints/brats_flair/"

    # 训练命令
    train_command = [
        "python", "pGAN.py",
        "--dataroot", dataroot,
        "--name", name,
        "--which_direction", which_direction,
        "--lambda_A", str(lambda_A),
        "--batchSize", str(batch_size),
        "--output_nc", str(output_nc),
        "--input_nc", str(input_nc),
        "--gpu_ids", str(gpu_ids),
        "--niter", str(niter),
        "--niter_decay", str(niter_decay),
        "--save_epoch_freq", str(save_epoch_freq),
        "--lambda_vgg", str(lambda_vgg),
        "--checkpoints_dir", checkpoints_dir,
        "--training"
    ]

    print("train...")
    subprocess.run(train_command)

def test():
    # 设置测试相关的参数
    dataroot = "datasets/REMBRANDT/"
    name = "pGAN_run"
    which_direction =  "AtoB"
    output_nc = 1
    input_nc = 3
    how_many = 50000
    results_dir = "result/REMBRANDT/"
    checkpoints_dir = "checkpoints/brats_t1ce/"

    # 测试命令
    test_command = [
        "python", "pGAN.py",
        "--dataroot", dataroot,
        "--name", name,
        "--which_direction", which_direction,
        "--phase", "test",
        "--output_nc", str(output_nc),
        "--input_nc", str(input_nc),
        "--how_many", str(how_many),
        "--results_dir", results_dir,
        "--checkpoints_dir", checkpoints_dir
    ]

    print("test...")
    subprocess.run(test_command)

if __name__ == "__main__":

    train()
    test()
    print("done")



