# Reciprocal Assistance Imputation-Representation Learning for Glioma Diagnosis with Incomplete MRI Sequences
A unified reciprocal assistance imputation-representation learning framework (REPAIR) is proposed for glioma (GM) diagnosis modeling with incomplete MRI sequences.


# REPAIR
The REPAIR adapts to varying scenarios with arbitrary missing of MRI sequences, rendering it a practical tool for handing comprehensive real-word medical data in complex clinical settings. REPAIR facilitates a cooperative process between missing value imputation and multi-sequence MRI fusion by leveraging existing samples to inform the imputation of missing values. This, in turn, facilitates the learning of a shared latent representation, which reciprocally guides more accurate imputation of missing values. To tailor the learned representation for downstream tasks, a novel ambiguity-aware intercorrelation regularization is introduced to equip REPAIR by correlating imputation ambiguity and its impacts conveying to the learned representation via a fuzzy paradigm. Additionally, a multimodal structural calibration constraint is devised to correct for the structural shift caused by missing data, ensuring structural consistency between the learned representations and the actual data. The proposed methodology is extensively validated on eight GM datasets with incomplete MRI sequences and six clinical datasets from other diseases with incomplete imaging modalities.
![image](https://github.com/user-attachments/assets/da41fd58-84da-44ac-b23a-3d927883e6ab)

# Model training & validation
One parameter was tuned at a time while keeping the others fixed. Varying missing rates ranging from 10% to 60% were simulated to validate the proposed method’s capability in handling incomplete MRI sequences. The learned representation was subsequently used to train a typical logistic regression classifier (with ‘l_2’ penalty while other parameters set as default) for all tasks across all evaluated datasets. We evaluated the method via a 5-fold cross-validation strategy.The classification performance was quantified by several metrics, including the area under curve (AUC), accuracy (ACC), balanced accuracy (BAcc), sensitivity (SEN) and specificity (SPE).

# Usage
Main file："main.m"

Input: Enter a training set (tr_X) and a testing set (tt_X) containing the missing sequences, and a label set (tr_Y) for the training set

Initialization: Data imputation matrix (Q), Latent representation matrix (V), Projection matrix (P), Reconstruction matrix (H)

Output: Latent representation matrix of training set (V_train) and Latent representation matrix of testing set (V_test)

# Note
Other investigators are invited to share their data in order to further improve the generalization capability of the current model. The model should only be used to support clinical diagnosis by health care professionals as a complementary tool for GM diagnosis with incomplete MRI sequences. Any responsibility for using this model and its results will rest solely by the health care professional using the model. Using it you should understand and agree that this tool is not responsible or liable for any claim, loss, or damage resulting from its use. While we try to keep the information on the tool as accurate as possible, we disclaim any warranty concerning its accuracy, timeliness, and completeness, and any other warranty, express or implied, including warranties of merchantability or fitness for a particular purpose.
