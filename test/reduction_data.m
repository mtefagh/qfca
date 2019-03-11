clear all; close all; clc;
bigg_universal = jsondecode(fileread('universal_model.json'));

reactions = bigg_universal.reactions;
id = strings(length(reactions), 1);
for i = 1:length(reactions)
    id(i) = reactions(i).id;
end

metabolites = bigg_universal.metabolites;
mets = strings(length(metabolites), 1);
for i = 1:length(metabolites)
    mets(i) = metabolites(i).id;
end 

S = sparse(length(metabolites), length(reactions));
for i = 1:length(reactions)
    metabolites = reactions(i).metabolites;
    mets_names = fieldnames(metabolites);
    for j = 1:length(mets_names)
        S(mets == mets_names{j}, i) = metabolites.(mets_names{j});
    end
end

model.S = S;
model.rev = zeros(length(id), 1);
model.rxns = cellstr(id);
model.mets = cellstr(mets);
clear i j S id mets mets_names metabolites reactions bigg_universal;