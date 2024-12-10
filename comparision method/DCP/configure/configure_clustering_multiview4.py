def get_default_config(data_name):
    if data_name in ['Caltech101-20']:
        return dict(
            type='CG',  # other: CV
            view=4,  # 设置为4视图
            seed=1,
            training=dict(
                lr=3.0e-4,
                start_dual_prediction=50,
                batch_size=256,
                epoch=100,
                alpha=10,
                lambda2=0.1,
                lambda1=0.1,
            ),
            Autoencoder=dict(
                arch1=[1984, 1024, 1024, 1024, 128],
                arch2=[512, 1024, 1024, 1024, 128],
                arch3=[928, 1024, 1024, 1024, 128],
                arch4=[1024, 1024, 1024, 1024, 128],  # 添加第四视图的Autoencoder结构
                activations='relu',
                batchnorm=True,
            ),
            Prediction=dict(
                arch1=[128, 256, 128],
                arch2=[128, 256, 128],
                arch3=[128, 256, 128],
                arch4=[128, 256, 128],  # 添加第四视图的预测结构
            ),
        )

    elif data_name in ['Scene_15']:
        """The default configs."""
        return dict(
            type='CG',
            seed=8,
            view=4,  # 设置为4视图
            training=dict(
                lr=1.0e-4,
                start_dual_prediction=100,
                batch_size=256,
                epoch=300,
                alpha=10,
                lambda2=0.1,
                lambda1=0.1
            ),
            Autoencoder=dict(
                arch1=[20, 1024, 1024, 1024, 128],
                arch2=[59, 1024, 1024, 1024, 128],
                arch3=[40, 1024, 1024, 1024, 128],
                arch4=[30, 1024, 1024, 1024, 128],  # 添加第四视图的Autoencoder结构
                activations='relu',
                batchnorm=True,
            ),
            Prediction=dict(
                arch1=[128, 256, 128],
                arch2=[128, 256, 128],
                arch3=[128, 256, 128],
                arch4=[128, 256, 128],  # 添加第四视图的预测结构
            ),
        )
    elif data_name in ['brats_hgg_lgg', 'brats_2', 'dataverse', 'REMBRANDT']:
        """The default configs."""
        return dict(
            type='CG',
            seed=4,
            view=4,  # 设置为4视图
            Autoencoder=dict(
                arch1=[109, 1024, 1024, 1024, 40],
                arch2=[109, 1024, 1024, 1024, 40],
                arch3=[109, 1024, 1024, 1024, 40],
                arch4=[109, 1024, 1024, 1024, 40],  # 添加第四视图的Autoencoder结构
                activations='relu',
                batchnorm=True,
            ),
            Prediction=dict(
                arch1=[128, 256, 128],
                arch2=[128, 256, 128],
                arch3=[128, 256, 128],
                arch4=[128, 256, 128],  # 添加第四视图的预测结构
            ),
            training=dict(
                lr=1.0e-4,
                start_dual_prediction=100,
                batch_size=256,
                epoch=500,
                alpha=10,
                lambda2=0.1,
                lambda1=0.1,
            ),
        )

    else:
        raise Exception('Undefined data name')
