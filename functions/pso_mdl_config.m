function cfg = pso_mdl_config(profile)
%PSO_MDL_CONFIG Define portable inputs and the validated research parameters.
%   CFG = PSO_MDL_CONFIG("formal") preserves the research search settings.
%   Edit the input/output paths below before running Steps 00-03.
%   The four original input groups are provided under data/private/.
if nargin < 1, profile = "formal"; end
profile = lower(string(profile));
assert(isscalar(profile) && profile == "formal", ...
    'The public configuration uses the "formal" profile.');
root = fileparts(fileparts(mfilename('fullpath')));
cfg = struct();
cfg.Profile = char(profile);
% Accept one descriptor TXT, or a directory of ordered TXT parts with
% descriptor_parts.json. Both represent the same complete candidate table.
cfg.Input.CandidateDescriptors = fullfile(root, 'data', 'private', 'candidate_descriptor_parts');
% cfg.Input.CandidateDescriptors = fullfile(root, 'data', 'private', 'candidates_descriptors.txt');
cfg.Input.ReferenceDescriptors = fullfile(root, 'data', 'private', 'reference_descriptors.txt');
cfg.Input.CandidateManifest = fullfile(root, 'data', 'private', 'candidate_manifest.tsv');
cfg.Input.ReferenceManifest = fullfile(root, 'data', 'private', 'reference_manifest.tsv');
cfg.Input.MissingnessCache = '';
cfg.Output.Directory = fullfile(root, 'results', 'research');
cfg.Output.Overwrite = false;

cfg.Identity.CandidateIdentifierColumn = "INCHIKEY";
cfg.Identity.ReferenceIdentifierColumn = "INCHIKEY";
cfg.Identity.TieBreakColumns = ["IDNUMBER", "DragonName"];
cfg.Identity.MissingIdentifierPolicy = "exclude";
cfg.Identity.ExcludeReferenceOverlap = true;
cfg.Background.SelectionMode = "property_matched";
cfg.Background.SizeMode = "capped_fraction";
cfg.Background.FixedCount = 20000;
cfg.Background.Fraction = 0.10;
cfg.Background.MaximumCount = 20000;
cfg.Background.MinimumCount = 100;
cfg.Background.MinimumHoldoutCount = 100;
cfg.Background.RandomSeed = 20260724;
cfg.Background.MatchingProperties = ["MW", "MLOGP", "TPSA(Tot)", "nHDon", "nHAcc", "RBN"];
cfg.Preprocessing.MaximumCandidateMissingRate = 0.05;
cfg.Preprocessing.WeightFactor = 0.2;
cfg.Preprocessing.NearZeroThreshold = 0.0001;
cfg.Preprocessing.ScaleSampleStride = 1;

% Historical NIRO runs used 28 computational threads and the AVX2 CNR branch.
% The branch is checked, not forced: MKL may retain AUTO on non-Intel CPUs.
cfg.Execution.ThreadLimit = 28;
cfg.Execution.ExpectedBLASBranch = "AVX2";

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
cfg.PSO.StallWindow = 200;
cfg.PSO.StallRelativeTolerance = 1e-4;
cfg.PSO.DistanceEpsilon = 1e-12;

cfg.Ranking.ExpectedReferenceIndex = [];
cfg.Ranking.StabilityTopK = 1000;
cfg.Ranking.StabilityCurveTopK = [20 50 100 200 500 1000 2000 5000 10000];
cfg.Ranking.FigureTopCandidateCount = 20;
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
cfg.IO.CandidateChunkRows = 2000;
cfg.IO.ReferenceChunkRows = 100;
cfg.IO.PreprocessingChunkRows = 500;
cfg.IO.ScoringChunkRows = 500;
cfg.IO.ProgressInterval = 25000;
cfg.IO.ScoringProgressInterval = 10000;
cfg.IO.ComputeInputSHA256 = true;

end
