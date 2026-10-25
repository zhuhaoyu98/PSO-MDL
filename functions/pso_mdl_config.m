function cfg = pso_mdl_config(profile)
%PSO_MDL_CONFIG Central settings for the five-step PSO-MDL workflow.
%   CFG = PSO_MDL_CONFIG returns the formal research configuration.
%   Edit this file for a new dataset or output directory; the Step files and
%   numerical helpers then use the same inputs and method parameters.
if nargin < 1, profile = "formal"; end
profile = lower(string(profile));
assert(isscalar(profile) && profile == "formal", ...
    'The public configuration uses the "formal" profile.');
root = fileparts(fileparts(mfilename('fullpath')));
cfg = struct();
cfg.Profile = char(profile);

%% 1. Four original input groups and one output directory
% Paths are relative to this project. The 100 TXT parts form one candidate table.
cfg.Input.CandidateDescriptors = fullfile(root, 'data', 'private', 'candidate_descriptor_parts');
% cfg.Input.CandidateDescriptors = fullfile(root, 'data', 'private', 'candidates_descriptors.txt');
cfg.Input.ReferenceDescriptors = fullfile(root, 'data', 'private', 'reference_descriptors.txt');
cfg.Input.CandidateManifest = fullfile(root, 'data', 'private', 'candidate_manifest.tsv');
cfg.Input.ReferenceManifest = fullfile(root, 'data', 'private', 'reference_manifest.tsv');
cfg.Input.MissingnessCache = '';
cfg.Output.Directory = fullfile(root, 'results', 'research');
cfg.Output.Overwrite = false;

%% 2. Compound identities and background selection
% Deduplicate by full INCHIKEY; this does not standardize molecular structures.
cfg.Identity.CandidateIdentifierColumn = "INCHIKEY";
cfg.Identity.ReferenceIdentifierColumn = "INCHIKEY";
cfg.Identity.TieBreakColumns = ["IDNUMBER", "DragonName"];
cfg.Identity.MissingIdentifierPolicy = "exclude";
cfg.Identity.ExcludeReferenceOverlap = true;
cfg.Background.SelectionMode = "property_matched";
% With capped_fraction: at most 20,000 and approximately 10% of eligible records.
% FixedCount is used only when SizeMode="fixed".
cfg.Background.SizeMode = "capped_fraction";
cfg.Background.FixedCount = 20000;
cfg.Background.Fraction = 0.10;
cfg.Background.MaximumCount = 20000;
cfg.Background.MinimumCount = 100;
cfg.Background.MinimumHoldoutCount = 100;
cfg.Background.RandomSeed = 20260724;
% RandomSeed applies to random mode; property matching is deterministic.
cfg.Background.MatchingProperties = ["MW", "MLOGP", "TPSA(Tot)", "nHDon", "nHAcc", "RBN"];

%% 3. Descriptor availability, imputation and scaling
% Medians and mean-AD divisors are fitted on background compounds only.
% WeightFactor multiplies selected divisors, not descriptor values or J.
cfg.Preprocessing.MaximumCandidateMissingRate = 0.05;
cfg.Preprocessing.WeightFactor = 0.2;
cfg.Preprocessing.NearZeroThreshold = 0.0001;
cfg.Preprocessing.ScaleSampleStride = 1;

%% 4. Numerical execution environment
% Historical NIRO runs used 28 computational threads and the AVX2 CNR branch.
% The branch is checked, not forced: MKL may retain AUTO on non-Intel CPUs.
cfg.Execution.ThreadLimit = 28;
cfg.Execution.ExpectedBLASBranch = "AVX2";

%% 5. Independent PSO searches of J = -R*F*P
% R=(D2-D1)/(D1+epsilon); F counts descriptors; P is a mean percentile (0-100).
cfg.PSO.Seeds = (0:63).';
cfg.PSO.PopulationSize = 50;
cfg.PSO.GenerationCount = 3000;
cfg.PSO.ScheduleHorizon = 3000;
cfg.PSO.InertiaWeightStart = 0.80;
cfg.PSO.InertiaWeightEnd = 0.35;
cfg.PSO.CognitiveCoefficient = 1.49618;
cfg.PSO.SocialCoefficient = 1.49618;
cfg.PSO.ParticleMutationRateStart = 0.10;
cfg.PSO.ParticleMutationRateEnd = 0.03;
cfg.PSO.ParticleMutationRateFloor = 0.003;
cfg.PSO.FeatureMutationRateStart = 0.008;
cfg.PSO.FeatureMutationRateEnd = 0.001;
cfg.PSO.FeatureMutationRateFloor = 8e-5;
cfg.PSO.MutationSchedulePower = 2.0;
cfg.PSO.MutationQuenchStartFraction = 0.85;
cfg.PSO.MutationNoiseScale = 0.10;
cfg.PSO.VelocityLimit = 4.0;
cfg.PSO.PositionLimit = 4.0;
% Stall settings summarize final-window change; they do not stop a run early.
cfg.PSO.StallWindow = 200;
cfg.PSO.StallRelativeTolerance = 1e-4;
cfg.PSO.DistanceEpsilon = 1e-12;

%% 6. Equal-weight rank aggregation and ranking-stability summaries
% Leave ExpectedReferenceIndex empty to use the fixed medoid determined from X_R.
cfg.Ranking.ExpectedReferenceIndex = [];
cfg.Ranking.StabilityTopK = 1000;
cfg.Ranking.StabilityCurveTopK = [20 50 100 200 500 1000 2000 5000 10000];
cfg.Ranking.FigureTopCandidateCount = 20;

%% 7. Figure 3/4 export and independent ensemble-size comparisons
cfg.Figures.Show = false;
cfg.Figures.Overwrite = false;
cfg.Figures.Dpi = 600;
cfg.Figures.FontName = 'Arial';
cfg.Figures.PublicationFontSize = 8.5;
cfg.Figures.FitnessDisplayScale = 1e6;
cfg.Figures.RegularSizeCm = [8.4 6.6];
cfg.Figures.TallSizeCm = [9.0 10.5];
cfg.Figures.EnsembleSubsetSizes = [4 8 16 32];
cfg.Figures.EnsembleRepeats = 100;
cfg.Figures.EnsembleRandomSeed = 20260808;
cfg.Figures.EnsembleTopK = 1000;

%% 8. Streaming settings (memory/progress only, not statistical sampling)
cfg.IO.CandidateChunkRows = 2000;
cfg.IO.ReferenceChunkRows = 100;
cfg.IO.PreprocessingChunkRows = 500;
cfg.IO.ScoringChunkRows = 500;
cfg.IO.ProgressInterval = 25000;
cfg.IO.ScoringProgressInterval = 10000;
cfg.IO.ComputeInputSHA256 = true;

end
