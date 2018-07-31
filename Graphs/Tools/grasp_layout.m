%Computes a layout of a graph using GraphViz.
%
%   layout = GRASP_LAYOUT(graph) given a graph, computes a layout
%   (i.e. a n-by-2 matrix) of the nodes.
%
%   GRASP_LAYOUT(..., options) extra options to control the process:
%
%   options.input_file: input file name
%
%   options.output_file: output file name
%
%   options.remove_files: boolean set to true if input_file and output_file
%       should be deleted at the end
%
%   options.algo: GraphViz algorithm (possible values are 'spring_basic'
%       (default), 'spring_force', 'spring_force_multiscale', 'tree',
%       'circular' and 'radial')
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Inspired by GraphViz4Matlab <https://code.google.com/p/graphviz4matlab/>

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
% 
% benjamin.girault@ens-lyon.fr
% benjamin.girault@usc.edu
% 
% This software is a computer program whose purpose is to provide a Matlab
% / Octave toolbox for handling and displaying graph signals.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.

function layout = grasp_layout(graph, varargin)
    %% Parameters
    default_param = struct(...
        'input_file', 'input.dot',...
        'output_file', 'output.dot',...
        'remove_files', true,...
        'algo', 'spring_basic');
    if nargin == 1
        options = struct;
    elseif nargin > 2
        options = cell2struct(varargin(2:2:end), varargin(1:2:end), 2);
    else
        options = varargin{1};
    end
    options = grasp_merge_structs(default_param, options);
    %% Initialisation
    is_directed = grasp_is_directed(graph);
    graph_type = 'graph';
    edge_text = '--';
    if is_directed
        graph_type = 'digraph';
        edge_text = '->';
    end
    
    %% Write dot file
    fid = fopen(options.input_file,'w');
    fprintf(fid, '%s G {\ncenter = 1;\nsize="10,10";\n', graph_type);
    N = grasp_nb_nodes(graph);
    for i = 1:N
        fprintf(fid,'%d;\n',i);
    end
    for i = 1:N
        begining = i + 1;
        if is_directed
            begining = 1;
        end
        for j = begining:N
            if i == j
                continue;
            end
            if graph.A(i,j) > 0
              fprintf(fid, '%d %s %d [weight=1, w=%f];\n', i, edge_text, j, full(graph.A(i, j)));
            end
        end 
    end
    fprintf(fid,'}');
    fclose(fid);
    
    %% Call GraphViz
    gviz_tool = 'neato';
    switch options.algo
        case 'tree'
            gviz_tool = 'dot';
        case 'circular'
            gviz_tool = 'circo';
        case 'radial'
            gviz_tool = 'twopi';
        case 'spring_force'
            gviz_tool = 'fdp';
        case 'spring_force_multiscale'
            gviz_tool = 'sfdp';
        case 'spring_basic'
        otherwise
            gviz_tool = 'neato';
    end
    err = system([gviz_tool, ' -Tdot -Gmaxiter=5000 -Gstart=7 -o ', options.output_file,' ', options.input_file]);
    if (err)
        if options.remove_files
            delete(options.input_file);
            delete(options.output_file);
        end
        error('Sorry, unknown GraphViz failure, please try another layout');
    end
    
    %% Parsing
    fid = fopen(options.output_file, 'r');
    % Read each block
    content = textscan(fid, '%s', 'delimiter', ';{}', 'EndOfLine', '');
    fclose(fid);
    content = content{:};
    % Remove newlines
    content = regexprep(content, '(\t|\n|\r)*', '');
    [~, bounds] = strtok(content{cellfun(@(x)~isempty(x),strfind(content,'graph [bb="'))},'"');
    bounds = textscan(strtok(bounds,'"'),'%n','delimiter',',');
    bounds = bounds{:}';
    content(cellfun(@(x)~isempty(x),strfind(content,sprintf(' %s ', edge_text))))=[]; % delete edge info, we don't need it
    content(cellfun(@(x)isempty(x),strfind(content,'pos')))=[];   % keep only positions
    regexprep(content, '^\d .*pos="%f,%f".*$', '$1,$2,$3');
    [~, remaining] = strtok(content,'"');
    [layout, ~] = strtok(remaining,'"');
    layout = cellfun(@(str)textscan(str,'%n','delimiter',','),layout,'UniformOutput',false);
    layout = [layout{:}];
    layout = [layout{:}]';
    
    %% Scaling in [0;10]²
    layout(:, 1) = (layout(:, 1) - bounds(1)) ./ (bounds(3) - bounds(1)) * 10;
    layout(:, 2) = (layout(:, 2) - bounds(2)) ./ (bounds(4) - bounds(2)) * 10;
    
    %% Clean up
    if options.remove_files
        delete(options.input_file);
        delete(options.output_file);
    end
end