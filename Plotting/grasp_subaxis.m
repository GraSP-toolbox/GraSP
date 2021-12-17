%Creates a sub axis within a figure.
%
%   ah = GRASP_SUBAXIS(nb_rows, nb_cols, cur_cell) creates the
%   (cur_cell)^th axes within an array of axes of size nb_rows x nb_cols.
%
%   ah = GRASP_SUBAXIS(fh, ...) uses the provided figure handle.
%
%   ah = GRASP_SUBAXIS(..., options)  optional parameters:
%
%   options.Spacing: array [SH SV]
%   options.SpacingXX:
%   options.SY: horizontal (XX = Horizontal / Y = H) or vertical
%       (XX = Vertical / Y = V) spacing between axes
%
%   options.Padding: array [PL PR PT PB]
%   options.PaddingXX:
%   options.PY: horizontal or [PL PR] (XX = Horizontal / Y = H),
%       or vertical or [PT PB] (XX = Vertical / Y = V), or left
%       (XX = Left / Y = L), or right (XX = Right / Y = R)
%       padding of an axe.
%
%   options.Margin: array [ML MR MT MB]
%   options.MarginXX:
%   options.MY: horizontal or [ML MR] (XX = Horizontal / Y = H),
%       or vertical or [MT MB] (XX = Vertical / Y = V), or left
%       (XX = Left / Y = L), or right (XX = Right / Y = R)
%       margin of the figure.
%
% Authors:
%  - Benjamin Girault <benjamin.girault@ens-lyon.fr>
%  - Benjamin Girault <benjamin.girault@usc.edu>
%  - Benjamin Girault <benjamin.girault@ensai.fr>
%  - Based on SubAxis (Similar API, re-implemented) (http://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot)

% Copyright Benjamin Girault, École Normale Supérieure de Lyon, FRANCE /
% Inria, FRANCE (2015)
% Copyright Benjamin Girault, University of Sourthern California, Los
% Angeles, California, USA (2016)
% Copyright Benjamin Girault, École Nationale de la Statistique et de
% l'Analyse de l'Information, Bruz, FRANCE (2020-2021)
% 
% benjamin.girault@ens-lyon.fr
% benjamin.girault@usc.edu
% benjamin.girault@ensai.fr
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

function ah = grasp_subaxis(varargin)
    %% Parameters
    param = struct(...
        'SpacingVertical', 0.05,...
        'SpacingHorizontal', 0.05,...
        'PaddingLeft', 0,...
        'PaddingRight', 0,...
        'PaddingTop', 0,...
        'PaddingBottom', 0,...
        'MarginLeft', .1,...
        'MarginRight', .1,...
        'MarginTop', .1,...
        'MarginBottom', .1);
    
    if nargin >= 4 && isnumeric(varargin{4})
        fh = varargin{1};
        nb_rows = varargin{2};
        nb_cols = varargin{3};
        cur_cell = varargin{4};
        if nargin >= 5
            if nargin == 5
                param = grasp_merge_structs(param, varargin{5}, false);
            else
                param = grasp_merge_structs(param, cell2struct(varargin(6:2:end), varargin(5:2:end), 2), false);
            end
        end
    else
        fh = gcf;
        nb_rows = varargin{1};
        nb_cols = varargin{2};
        cur_cell = varargin{3};
        if nargin >= 4
            if nargin == 4
                param = grasp_merge_structs(param, varargin{4}, false);
            else
                param = grasp_merge_structs(param, cell2struct(varargin(5:2:end), varargin(4:2:end), 2), false);
            end
        end
    end
    
    if isfield(param, 'Spacing')
        param.S = param.Spacing;
    end
    if isfield(param, 'S')
        if numel(param.S) == 1
            param.SpacingHorizontal = param.S;
            param.SpacingVertical = param.S;
        else
            param.SpacingHorizontal = param.S(1);
            param.SpacingVertical = param.S(2);
        end
    end
    if isfield(param, 'SV')
        param.SpacingVertical = param.SV;
    end
    if isfield(param, 'SH')
        param.SpacingHorizontal = param.SH;
    end
    
    if isfield(param, 'Padding')
        param.P = param.Padding;
    end
    if isfield(param, 'P')
        if numel(param.P) == 1
            param.PaddingLeft = param.P;
            param.PaddingRight = param.P;
            param.PaddingTop = param.P;
            param.PaddingBottom = param.P;
        elseif numel(param.P) == 2
            param.PaddingLeft = param.P(1);
            param.PaddingRight = param.P(1);
            param.PaddingTop = param.P(2);
            param.PaddingBottom = param.P(2);
        else
            param.PaddingLeft = param.P(1);
            param.PaddingRight = param.P(2);
            param.PaddingTop = param.P(3);
            param.PaddingBottom = param.P(4);
        end
    end
    if isfield(param, 'PaddingHorizontal')
        param.PH = param.PaddingHorizontal;
    end
    if isfield(param, 'PH')
        if numel(param.PH) == 1
            param.PaddingLeft = param.PH;
            param.PaddingRight = param.PH;
        else
            param.PaddingLeft = param.PH(1);
            param.PaddingRight = param.PH(2);
        end
    end
    if isfield(param, 'PaddingVertical')
        param.PV = param.PaddingVertical;
    end
    if isfield(param, 'PV')
        if numel(param.PV) == 1
            param.PaddingTop = param.PV;
            param.PaddingBottom = param.PV;
        else
            param.PaddingTop = param.PV(1);
            param.PaddingBottom = param.PV(2);
        end
    end
    if isfield(param, 'PL')
        param.PaddingLeft = param.PL;
    end
    if isfield(param, 'PR')
        param.PaddingRight = param.PR;
    end
    if isfield(param, 'PT')
        param.PaddingTop = param.PT;
    end
    if isfield(param, 'PB')
        param.PaddingBottom = param.PB;
    end
    
    if isfield(param, 'Margin')
        param.M = param.Margin;
    end
    if isfield(param, 'M')
        if numel(param.M) == 1
            param.MarginLeft = param.M;
            param.MarginRight = param.M;
            param.MarginTop = param.M;
            param.MarginBottom = param.M;
        elseif numel(param.M) == 2
            param.MarginLeft = param.M(1);
            param.MarginRight = param.M(1);
            param.MarginTop = param.M(2);
            param.MarginBottom = param.M(2);
        else
            param.MarginLeft = param.M(1);
            param.MarginRight = param.M(2);
            param.MarginTop = param.M(3);
            param.MarginBottom = param.M(4);
        end
    end
    if isfield(param, 'MarginHorizontal')
        param.MH = param.MarginHorizontal;
    end
    if isfield(param, 'MH')
        if numel(param.MH) == 1
            param.MarginLeft = param.MH;
            param.MarginRight = param.MH;
        else
            param.MarginLeft = param.MH(1);
            param.MarginRight = param.MH(2);
        end
    end
    if isfield(param, 'MarginVertical')
        param.MV = param.MarginVertical;
    end
    if isfield(param, 'MV')
        if numel(param.MV) == 1
            param.MarginTop = param.MV;
            param.MarginBottom = param.MV;
        else
            param.MarginTop = param.MV(1);
            param.MarginBottom = param.MV(2);
        end
    end
    if isfield(param, 'ML')
        param.MarginLeft = param.ML;
    end
    if isfield(param, 'MR')
        param.MarginRight = param.MR;
    end
    if isfield(param, 'MT')
        param.MarginTop = param.MT;
    end
    if isfield(param, 'MB')
        param.MarginBottom = param.MB;
    end
    
    %% Cell indices
    i = floor((cur_cell - 1) / nb_cols) + 1;
    j = mod((cur_cell - 1), nb_cols) + 1;
    
    %% Axes boundaries
    SH = param.SpacingHorizontal;
    SV = param.SpacingVertical;
    MH = param.MarginLeft + param.MarginRight;
    MV = param.MarginTop + param.MarginBottom;
    cellwidth = ((1 - MH) - (nb_cols - 1) * SH) / nb_cols;
    cellheight = ((1 - MV) - (nb_rows - 1) * SV) / nb_rows;
    
    ML = param.MarginLeft;
%     MR = param.MarginRight;
%     MT = param.MarginTop;
    MB = param.MarginBottom;
    PL = param.PaddingLeft;
    PR = param.PaddingRight;
    PT = param.PaddingTop;
    PB = param.PaddingBottom;
    xleft = ML + cellwidth * (j - 1) + SH * (j - 1) + PL;
    xwidth = cellwidth - PL - PR;
    ybottom = MB + cellheight * (nb_rows - i) + SV * (nb_rows - i) + PB;
    yheight = cellheight - PB - PT;

    %% Axes
    figure(fh);
    ah = subplot('Position', [xleft, ybottom, xwidth, yheight]);
end
