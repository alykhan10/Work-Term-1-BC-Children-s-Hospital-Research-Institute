%% Defining the Matrix Tables.
ZTABLE = readtable("C:\Users\alykh\Downloads\Woodward Lab\Constrained Principal Component Analysis (CPCA) Class\SPQ_ABQ_BFAS_for_CPCA.csv");
GTABLE = readtable("C:\Users\alykh\Downloads\Woodward Lab\Constrained Principal Component Analysis (CPCA) Class\SPQ_ABQ_BFAS_for_CPCA.csv");
ZTABLE = removevars(ZTABLE, ["BFAS_2","BFAS_3","BFAS_4","BFAS_5","BFAS_6","BFAS_7","BFAS_8","BFAS_9","BFAS_10","BFAS_11","BFAS_12","BFAS_13","BFAS_14","BFAS_15","BFAS_16","BFAS_17","BFAS_18","BFAS_19","BFAS_20","BFAS_21","BFAS_22","BFAS_23","BFAS_24","BFAS_25","BFAS_26","BFAS_27","BFAS_28","BFAS_29","BFAS_30","BFAS_31","BFAS_32","BFAS_33","BFAS_34","BFAS_35","BFAS_36","BFAS_37","BFAS_38","BFAS_39","BFAS_40","BFAS_41","BFAS_42","BFAS_43","BFAS_44","BFAS_45","BFAS_46","BFAS_47","BFAS_48","BFAS_49","BFAS_50","BFAS_51","BFAS_52","BFAS_53","BFAS_54","BFAS_55","BFAS_56","BFAS_57","BFAS_58","BFAS_59","BFAS_60","BFAS_61","BFAS_62","BFAS_63","BFAS_64","BFAS_65","BFAS_66","BFAS_67","BFAS_68","BFAS_69","BFAS_70","BFAS_71","BFAS_72","BFAS_73","BFAS_74","BFAS_75","BFAS_76","BFAS_77","BFAS_78","BFAS_79","BFAS_80","BFAS_81","BFAS_82","BFAS_83","BFAS_84","BFAS_85","BFAS_86","BFAS_87","BFAS_88","BFAS_89","BFAS_90","BFAS_91","BFAS_92","BFAS_93","BFAS_94","BFAS_95","BFAS_96","BFAS_97","BFAS_98","BFAS_99","BFAS_100"]);
ZTABLE = removevars(ZTABLE, ["ABQ_response_2","ABQ_response_3","ABQ_response_4","ABQ_response_5","ABQ_response_6","ABQ_response_7","ABQ_response_8","ABQ_response_9","ABQ_response_10","ABQ_response_11","ABQ_response_12","ABQ_response_13","ABQ_response_14","ABQ_response_15","ABQ_response_16","ABQ_response_17","ABQ_response_18","ABQ_response_19","ABQ_response_20","ABQ_response_21","ABQ_response_22","ABQ_response_23","ABQ_response_24","ABQ_response_25","ABQ_response_26","ABQ_response_27","ABQ_response_28","ABQ_response_29","ABQ_response_30","ABQ_response_31","ABQ_response_32","ABQ_response_33","ABQ_response_34","ABQ_response_35","ABQ_response_36","ABQ_response_37","ABQ_response_38","ABQ_response_39","ABQ_response_40"]);
GTABLE = removevars(GTABLE, ["SPQ_2","SPQ_3","SPQ_4","SPQ_5","SPQ_6","SPQ_7","SPQ_8","SPQ_9","SPQ_10","SPQ_11","SPQ_12","SPQ_13","SPQ_14","SPQ_15","SPQ_16","SPQ_17","SPQ_18","SPQ_19","SPQ_20","SPQ_21","SPQ_22","SPQ_23","SPQ_24","SPQ_25","SPQ_26","SPQ_27","SPQ_28","SPQ_29","SPQ_30","SPQ_31","SPQ_32","SPQ_33","SPQ_34","SPQ_35","SPQ_36","SPQ_37","SPQ_38","SPQ_39","SPQ_40","SPQ_41","SPQ_42","SPQ_43","SPQ_44","SPQ_45","SPQ_46","SPQ_47","SPQ_48","SPQ_49","SPQ_50","SPQ_51","SPQ_52","SPQ_53","SPQ_54","SPQ_55","SPQ_56","SPQ_57","SPQ_58","SPQ_59","SPQ_60","SPQ_61","SPQ_62","SPQ_63","SPQ_64","SPQ_65","SPQ_66","SPQ_67","SPQ_68","SPQ_69","SPQ_70","SPQ_71","SPQ_72","SPQ_73","SPQ_74","ABQ_response_2","ABQ_response_3","ABQ_response_4","ABQ_response_5","ABQ_response_6","ABQ_response_7","ABQ_response_8","ABQ_response_9","ABQ_response_10","ABQ_response_11","ABQ_response_12","ABQ_response_13","ABQ_response_14","ABQ_response_15","ABQ_response_16","ABQ_response_17","ABQ_response_18","ABQ_response_19","ABQ_response_20","ABQ_response_21","ABQ_response_22","ABQ_response_23","ABQ_response_24","ABQ_response_25","ABQ_response_26","ABQ_response_27","ABQ_response_28","ABQ_response_29","ABQ_response_30","ABQ_response_31","ABQ_response_32","ABQ_response_33","ABQ_response_34","ABQ_response_35","ABQ_response_36","ABQ_response_37","ABQ_response_38","ABQ_response_39","ABQ_response_40"]);

%% Obtaining the Matrix Array from the Table.
Z = table2array(ZTABLE(:,2:74));
G = table2array(GTABLE(:,2:100));

%% Performing Multivariate Multiple Regression.
Z_STANDARDIZED = (Z-(ones(size(Z)).*mean(Z)))./(ones(size(Z)).*std(Z));
C = inv(G'*G)*G'*Z_STANDARDIZED;
Z_HAT = G * C;
E = Z_STANDARDIZED - Z_HAT;

%% Performing Principal Component Analysis.
% Reducing the Matrices Using Singular Value Decomposition.
[U_Z,D_Z,V_Z] = svd(Z_STANDARDIZED,"econ");
[U_Z_HAT,D_Z_HAT,V_Z_HAT] = svd(G*C,"econ");
[U_E,D_E,V_E] = svd(E,"econ");

% Obtaining Eigenvalues from Singular Value Decomposition Matrix
EIGS_Z = (diag(D_Z).^2)./(length(Z_STANDARDIZED)-1);
EIGS_Z_HAT = (diag(D_Z_HAT).^2)./(length(Z_STANDARDIZED)-1);
EIGS_E = (diag(D_E).^2)./(length(Z_STANDARDIZED)-1);

U_Z = U_Z(:,1:2);
D_Z = D_Z(1:2,1:2);
V_Z = V_Z(:,1:2);

U_Z_HAT = U_Z_HAT(:,1:2);
D_Z_HAT = D_Z_HAT(1:2,1:2);
V_Z_HAT = V_Z_HAT(:,1:2);

U_E = U_E(:,1:2);
D_E = D_E(1:2,1:2);
V_E = V_E(:,1:2);

% Obtaining the Component Loadings for the Matrices.
COMPONENT_LOADINGS_Z = (V_Z*D_Z)./sqrt(length(Z_STANDARDIZED));
COMPONENT_LOADINGS_Z_HAT = (V_Z_HAT*D_Z_HAT)./sqrt(length(Z_HAT));
COMPONENT_LOADINGS_E = (V_E*D_E)./sqrt(length(E));

% Obtaining the Predictor Loadings for the Matrices.
PREDICTOR_LOADINGS = zeros(width(G),2);

for i = 1:width(G)
 for ii = 1:2
    PREDICTOR = corrcoef(G(:,i),U_Z_HAT(:,ii));
    predictor = PREDICTOR(1,2);

    PREDICTOR_LOADINGS(i,ii) = predictor;
 end
end

% Creating the CPCA Table.
[CPCA_Variance_table,CPCA_Percent_of_Source_table,CPCA_Percent_of_Overall_table] = Create_CPCA_Tables(Z_STANDARDIZED,G,C,E,2,EIGS_Z,EIGS_Z_HAT,EIGS_E);

%% Optimizing for Variance
% Performing the Varimax Rotation
[Component_Loadings_Overall_Rotated,Component_Loadings_Overall_Rotated_T] = varim(COMPONENT_LOADINGS_Z);
[Component_Loadings_Predicted_Rotated,Component_Loadings_Predicted_Rotated_T] = varim(COMPONENT_LOADINGS_Z_HAT);
[Component_Loadings_Residual_Rotated,Component_Loadings_Residual_Rotated_T] = varim(COMPONENT_LOADINGS_E);
[Component_Loadings_Predictor_Rotated,Component_Loadings_Predictor_Rotated_T] = varim(PREDICTOR_LOADINGS);

% Plotting the Rotated Component & Predictor Loadings
[Overall_Plot_Rotated,Predicted_Plot_Rotated,Residual_Plot_Rotated] = Component_Loading_2D_Plots(Variable_Labels,Predictor_Labels,Component_Loadings_Overall_Rotated,Component_Loadings_Predicted_Rotated,Component_Loadings_Residual_Rotated,Component_Loadings_Predictor_Rotated);