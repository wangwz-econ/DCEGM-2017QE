% This Matlab class cls_polyline implements a set of tools working with linearly interpolated functions.

% The original codes come from the replication code of paper, "The Endogenous Grid Method for 
%   Discrete-Continuous Dynamic Choice Models with (or without) Taste Shocks" (QE, 2017).

% Some details are taylored to fit my own programming habits and specific usage in script files
%   model2_DCEGM_deterministic_retirement.m.

classdef cls_polyline
  % This class defines a function linearly interpolated on a grid and implements the tools for 
  %   working with arrays of these lines.

  % Objects of this class can be put in arrays and manipulated in bulk, which means that each cls_polyline
  %   object in the array will conduct the following methods elementwise.

  % Methods: len() to return the number of points
  %					 sort() to sort the points
  %					 interpolate() to compute the function
  %					 upper_envelope() to find upper envelope of a set of lines
  %					 insert() to insert a new point
  %					 chop() to break one line into two
  % 				 grow() to join multiple lines together
  %					 thinout() to delete points by index
  %					 diff() to compare two polylines
  %					 plot() to make a plot of the function

  properties (Access=public)
    x = [];    % x-grid points
    y = [];    % function values
    label = "" % label of this line
  end

  methods
    %----------------------------------------------
    function obj = cls_polyline(varargin) 
      % Initialize the object
      if nargin>0
        obj.x = reshape(varargin{1}, 1, []);
      end
      if nargin>1
        obj.y = reshape(varargin{2}, 1, []);
      end
      if nargin>2 && ischar(varargin{3})
        obj.label = varargin{3};
      end
    end
    %----------------------------------------------
    function res = len(obj)
      % Number of grid points in this line/function.
      % When input obj is a cls_polylin array, it returns an array with the same size with each entry
      %   stores the length of x and y in the same entry of obj array.
      res = NaN(1, numel(obj));
      for k = 1:numel(obj)
        res(k) = numel(obj(k).x);
        if res(k)~=numel(obj(k).y)
          error('mismatch of the number of points: x has %d while y has %d elements',numel(obj(k).x),numel(obj(k).y))
        end
      end % end the loop for all elements in the obj array
      res = reshape(res, size(obj));
    end
    %----------------------------------------------
    function res = sort(obj)
      % This function sorts the polylines on x grid
      %   returns the sorted polylines, original obj is not modified
      res = NaN(1, numel(obj));
      for i = 1:numel(obj)
        [~, i1] = sort(obj(i).x); % obj(i).x is a vector
        res(i) = cls_polyline(obj(i).x(i1), obj(i).y(i1), obj(i).label);
      end % end the loop for all elements in the obj array
      res = reshape(res, size(obj));
    end
    %----------------------------------------------
    function [res, extrapflag] = interpolate(obj, vx)
      % This function computes the value in given points vx
      %   linear extrapolation allowed, indicator is created if nargout==2.
      % This function can take an object array as an input, 
      %   but vx (which contains the query points) must be a vector.
      % If an object input is passed to this function, then it returns an array, with rows being the 
      %  the original object array (the ith row is the ith linear index) and columns being the interpolated 
      %  values for vx (the jth col is the corresponding interpolated value for vx(j).
      if numel(obj)==1
        try
          res = interp1(obj.x, obj.y, vx, "linear", "extrap");
        catch ME
          steps = reshape(obj.x(2:end)-obj.x(1:end-1),[],1);
          fprintf('\nInterpolation error!\npolyline dims are (%d,%d)\nRepeated points: %d\nInverted steps: %d\n', ...
            size(obj, 1), size(obj, 2), sum(abs(steps))<eps, sum(steps<0));
          fprintf('Identifier: %s\nCall stack: ', ME.identifier);
          for i = 0:(numel(ME.stack)-1)
            if i>0
              fprintf(' >>> ')
            end
            fprintf('%s (line %d)', ME.stack(end-i).name, ME.stack(end-i).line);
          end
          fprintf('\n')
          fprintf('Interactive Keyboard Mode\nVariable steps contains the steps of the grid\n');
          keyboard
        end % end the try-catch statement
        if nargout >  1
          extrapflag = ~(vx>=min(obj.x) & vx<=max(obj.x)); % returns true if extrap is needed
        end
      else
        % interpolate each polyline if there are many output matrices with rows for each polyline,
        % columns of interpolated values
        xx = reshape(vx, 1, []);          % turn vx into a row vector
        res = NaN(numel(obj), numel(xx));
        for k=1:numel(obj)
          try
            res(k, 1:numel(xx)) = interp1(obj(k).x, obj(k).y, xx, 'linear', 'extrap');
          catch ME
            steps=@(nn) reshape(obj(nn).x(2:end)-obj(nn).x(1:end-1),[],1);
            fprintf('\n\nInterpolation error!\nk=%d, polyline dims are (%d,%d)\nRepeated points: %d\nInverted steps: %d\n', ...
              k,size(obj,1),size(obj,2),sum(abs(steps(k))<eps),sum(steps(k)<0));
            fprintf('Identifier: %s\nCall stack: ',ME.identifier);
            for i=0:numel(ME.stack)-1
              if i>0
                fprintf(' >>> ');
              end
              fprintf('%s (line %d)',ME.stack(end-i).name,ME.stack(end-i).line)
            end
            fprintf('\n')
            fprintf('INTERACTIVE KEYBOARD MODE\nUse steps(i) function to analyze the grids\n');
            keyboard
          end
          if nargout>1
            extrapflag = NaN(numel(obj), numel(xx));
            extrapflag(k,1:numel(xx)) = ~( xx>=min(obj(k).x) & xx<=max(obj(k).x) );
          end % end the if statement to check nargout
        end % end the for loop to go through all elements in the object array
      end % end the if-else statement to check the number of elements in the obj input
    end % end the function
    %----------------------------------------------
    function res=grow(obj, in2, front)
      % This function adds polylin2 to the end (or the front) of ojb.
      % Many-to-many extending is allowed.
      % Both obj and in2 can be cls_polyline arrays, but they must have the same size. And in this 
      %   case, the function is element-wise for each corresponding entry of obj and in2.
      if nargin<3
        front = false;
      end
      res = NaN(size(obj));
      for k1 = 1:numel(obj)
        if front
          res(k1) = cls_polyline( ...
            [reshape(in2(k1).x, 1, []), reshape(obj(k1).x, 1, [])], ...
            [reshape(in2(k1).y, 1, []), reshape(obj(k1).y, 1, [])], ...
            sprintf('%s + %s', in2(k1).label, obj(k1).label));
        else
          res(k1) = cls_polyline( ...
            [reshape(obj(k1).x,1,[]), reshape(in2(k1).x, 1, [])], ...
            [reshape(obj(k1).y,1,[]), reshape(in2(k1).y, 1, [])], ...
            sprintf('%s + %s',obj(k1).label,in2(k1).label));
        end % end the if statement to check front argument
      end % end the loop to go through all elements in the obj array
    end % end the function
    %----------------------------------------------
    function res = insert(obj, x, y, j)
      % This function inserts a point (x, y) after position j in the current grid, and j=0 adds the 
      %   point (x, y) to the front.
      % If x and y are vectors, multiple points are inserted.
      % If j is missing, the point is added in the end.
      % Inset into each polyline if there are many, i.e., obj is an cls_polyline array.
      % x, y can only be vectors or scalars.
      % Cases in the function to minimize overall runtime
      res = obj;
      if numel(x)>1 % here x and y are two vectors
        mx = numel(x);
        my = numel(y);
        if mx ~= my
          warning('The x and y inputs do not have the same length!')
        end
        for k = 1:numel(obj)
          if nargin<=3 || j>=numel(obj(k).x)
            res(k).x(end+1:end+mx) = x;
            res(k).y(end+1:end+my) = y;
          else
            res(k).x(j+mx+1:end+mx) = obj(k).x(j+1:end);
            res(k).x(j+1:j+mx) = x(:);
            res(k).y(j+my+1:end+my)=obj(k).y(j+1:end);
            res(k).y(j+1:j+my) = y(:);
          end % end the if-else statement to check the inserting position index j
        end % end the loop to go through all elements in the obj array
      elseif numel(obj)>1 % here the (x,y) are two scalars, obj is a cls_polyline object array
        for k=1:numel(obj)
          if nargin<=3 || j>=numel(obj(k).x)
            res(k).x(end+1) = x;
            res(k).y(end+1) = y;
          else
            res(k).x(j+2:end+1) = obj(k).x(j+1:end);
            res(k).x(j+1) = x;
            res(k).y(j+2:end+1) = obj(k).y(j+1:end);
            res(k).y(j+1) = y;
          end % end the if-else statement to check the inserting position index j
        end % end the loop to go through all elements in the obj array
      else % here the (x,y) are two scalars, obj is a cls_polyline object scalar
        if nargin<=3 || j>=numel(obj.x)
          res.x(end+1)=x;
          res.y(end+1)=y;
        else
          res.x(j+2:end+1)=obj.x(j+1:end);
          res.x(j+1)=x;
          res.y(j+2:end+1)=obj.y(j+1:end);
          res.y(j+1)=y;
        end % end the if-else statement to check the inserting position index j
      end % end the overall if-else statement to check whethter obj is an array, x and y are vectors
    end
    %----------------------------------------------
    function res = thinout(obj, index)
      % This function removes the indexed points from obj.
      res = obj;
      for k = 1:numel(res)
        ii = intersect(1:numel(res(k).x), index);
        res(k).x(ii) = [];
        res(k).y(ii) = [];
      end    
    end
    %----------------------------------------------
    function indx = diff(obj, pl2, significance)
      % This function returns the indexes of points in obj that are not in pl2.
      % obj can have multiple elements (obj can be an array), while pl2 is treated as scalar cls_polyline.
      if nargin < 3
        significance = 5; % equality is measured up to 10^-signif
      end
      x1 = round(pl2(1).x*(10^significance)) * 10^(-significance);
      y1 = round(pl2(1).y*(10^significance)) * 10^(-significance);
      for k = 1:numel(obj)
        a = round(obj(k).x*(10^significance)) * 10^(-significance);
        b = round(obj(k).y*(10^significance)) * 10^(-significance);
        indx{k} = find(~ismember(a, x1) | ~ismember(b, y1));
      end
      if numel(indx)==1
        indx = indx{1};
      else
        indx = reshape(indx, size(obj));
      end
    end
    %----------------------------------------------
    function [res1, res2] = chop(obj, j, repeat)
      % This function separates the grid into 1,..., j and j+1,..., end two parts.
      % If repeat=true, the boundary points are repeated in both resulting polylines
      if nargin<2
        error('Have to have one input for .chop')
      end
      for k=1:numel(obj) % chop each polyline if there are many
        if j>obj(k).len
          warning('Producing empty polyline by chopping at index j>len')
          j = obj(k).len;
        end
        res1(k) = cls_polyline(obj(k).x(1:j), obj(k).y(1:j), sprintf('%s (1:%d)', obj(k).label, j));
        if nargout > 1
          if repeat
            res2(k) = cls_polyline(obj(k).x(j:end), obj(k).y(j:end), ...
                                sprintf('%s (%d:%d)', obj(k).label, j, obj(k).len));
          else
            res2(k) = cls_polyline(obj(k).x(j+1:end), obj(k).y(j+1:end), ...
                                sprintf('%s (%d:%d)', obj(k).label, j+1, obj(k).len));
          end % end the if-else statement to check whether repeat==true
        end % end the if-else statement to check the number of output arguments
      end % end the for loop to go through all elements in obj array
      res1 = reshape(res1, size(obj));
      if nargout > 1
        res2 = reshape(res2, size(obj));
      end
    end
    %----------------------------------------------
    function [res, intersections] = upper_envelope(obj, fullinterval)
      % This function computes the upper envelope over the array of polylines.
      % It assumes that all grids are sorted or should be sorted, and treats them as sorted.
      % By default, only the overlapping segments are used for upper envelope calculation.
      % When fullinterval=true, all polylines are extended to union of invervals.
      if numel(obj)==1
        warning('Upper envelope is meant for an array of polylines')
        res = obj;
        if nargout>1
          intersections = cls_polyline([], [], "intersection points");
        end
        return
      end % end the if statement to check the number of elements in obj array

      l = obj.len;                % check that x and y are of the same size
      obj = obj(l>0);             % disregard all polylines of zero length
      pt = sort(unique([obj.x])); % collect all the x points in sorted row vector
      % IMPORTANT: When obj is a cls_polyline array, ojb.x returns all obj(k).x as vectors

      % interpolate all lines on all points recording the extrapolation cases
      [intr, extrflag] = obj.interpolate(pt);
      if ~fullinterval
        % disregard points where at least one line is extrapolated
        mask = sum(extrflag, 1) > 0;
        intr(:, mask) = [];
        pt(mask) = [];
        n = sum(~mask); % number of points in the overlap region
      else
        % disregard only points where particular lines are extrapolated (full interval!)
        intr(extrflag) = -Inf;
        n = numel(pt);
      end

      % find lines on the top
      maxintr = repmat(max(intr), size(intr,1), 1); % It is an array with each row being the original obj elements and columns being the maximized value when interpolate the union the obj(k).x
      top = (intr == maxintr);

      % build up the upper envelope
      res = cls_polyline(pt(1), maxintr(1,1), 'upper envelope');
      if nargout>1
        intersections = cls_polyline([], [], "intersection points");
      end
      k0 = find(top(:,1), 1, 'first'); % index of top line
      
      % loop through all points
      for j = 2:n
        k1 = find(top(:,j), 1, "first"); % index of next top line
        if k1 ~= k0
          % check if there is an intersection point between the lines:
          ln1=k0; ln2=k1; % intersections between these points
          pt1 = pt(j-1); pt2 = pt(j); % which lies between these points
          [y1, extr1] = obj(ln1).interpolate([pt1 pt2]);
          [y2, extr2] = obj(ln2).interpolate([pt1 pt2]); % and these function values (maybe extrapolated)
          % check that neither is extrapolated in both points
          % and that intersection point is inside the inverval <= function values are different at the borders
          if all(~[extr1 extr2]) && all(abs(y1-y2)>0)
            % find the intersection point or points
            while true
              pt3 = fzero( @(x) obj(ln2).interpolate(x) - obj(ln1).interpolate(x), [pt1 pt2] );
              pt3f = obj(ln1).interpolate(pt3);
              % check if there are lines above the found intersection
              [intr2, exrt2] = obj.interpolate(pt3); % interpolate all lines in the new point
              intr2(exrt2) = -Inf; % disregard the extrapolated points
              maxintr2 = repmat(max(intr2), size(intr2, 1, 1));
              ln3 = find(intr2==maxintr2, 1, 'first');
              if ln3 == ln1 || ln3==ln2
                % there are no other functions above!
                % add the intersection point
                res = res.inset(pt3, pt3f, res.len);
                if narargout > 1
                  intersections = intersections.inset(pt3, pt3f);
                end
                % maybe there are some more intersections before next point?
                if ln2 == k1
                  % no actually, because the left line is the one we started with
                  break
                else
                  % indeed, so update the interval of new search
                  ln1 = ln2;
                  pt1 = pt3;
                  ln2 = k1;
                  pt2 = pt(j);
                end
              else
                % there is line ln3 above the found intersection point
                % so, it is not on the upper envelope
                % need to search again
                ln2 = ln3; % new candidate
                pt2 = pt3;
              end
            end
          end
        end
        if any(abs(obj(k1).x-pt(j)) < eps) || j==n
          % add new grid point to the end
          res = res.inset(pt(j), maxintr(1, j));
        end
        k0 = k1; % next step
      end
    end

    function [res, indxremoved, newdots]=secondary_envelope(obj)
  		% This function computes the secondary envelope of the polyline.
  		% It returns cleaned polylines, indexes of removed points and new points as a polyline.
  		% In case of many polylines, clean up each.
  		% indxremoved is cell array when obj has many elements.
  		for k=1:numel(obj)
  			cur=obj(k);%current line
  			%identify the loop-back regions
  			ii=cur.x(2:end)>cur.x(1:end-1); %zeros for the loopback regions
  			sect=polyline;%sections
  			i=1;
  			while true
  				j=find(ii~=ii(1),1,'first');
  				if isempty(j)
  					%exit the loop if there are no more loop-backs
  					if i>1
  						%if sections already started, add the last one
  						sect(i)=cur;
  					end
  					break;
  				end
  				[sect(i), cur]=cur.chop(j,true); %true to repeat the boundary points
  				ii(1:j-1)=[];%chop ii array accordingly
  				i=i+1;
  			end
  			% perform secondary envelope if sections created
  			if numel(sect)>1
  				sect=sect.sort; %sort all sections since half of them are in opposite direction
  				[res(k), newdots(k)]=sect.upper_envelope(true); %true for full interval
  				%removed points indexes
  				indxremoved{k}=obj(k).diff(res(k),10); %index of dots in obj(k) but not in res(k)
  			else
  				%without loopbacks -- just copy the input
  				res(k)=obj(k);
  				indxremoved{k}=[];
  				newdots(k)=polyline;
  			end
  		end %next k
  		%return in same dimensions
  		res=reshape(res,size(obj));
  		indxremoved=reshape(indxremoved,size(obj));
  		if numel(obj)==1
  			%for a single polyline return indx as a vector
  			indxremoved=indxremoved{1};
  		end
  	end
  	%----------------------------------------------
  	function res=plot(obj,ax,varargin)
  		% This function plots the function on the given axes with optional line properties.
  		% It returns line handle(s)
  		% It allows additional input arguments as line properties.
  		% Additional arguments may also contain 'sort' to sort the data
  		l=obj.len;%check that x and y are of the same size
  		if all(l==0)
  			warning 'Nothing to plot'
  			return
  		end
  		if nargin>1
  			if isempty(ax)
  				fig1=figure('Color','white');
  				ax=axes('Parent',fig1);
  			elseif ~ishandle(ax)
  				error 'First argument should be axes handle or []'
  			end
  		end
  		if nargin>2
  			mask=cellfun(@(x) ischar(x)&strcmp(x,'sort'),varargin);%find 'sort'
  			if sum(mask)>0
  				dosort=true;
  			else
  				dosort=false;
  			end
  			opts=varargin(~mask);
  		else
  			opts={'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','MarkerEdgeColor','auto'};
  			dosort=false;%no sorting by default
  		end
  		for k=1:numel(obj) %draw all polylines if there are many
  			if l(k)>0
  				if dosort %sorting the points of the grid
  					[xx, j]=sort(obj(k).x);
  					xx=reshape(xx,[],1);
  					yy=reshape(obj(k).y(j),[],1);
  				else
  					xx=reshape(obj(k).x,[],1);
  					yy=reshape(obj(k).y,[],1);
  				end
  				if exist('ax')
  					hold(ax,'all');
  					res(k)=plot(ax,xx,yy,opts{:},'DisplayName',obj(k).label);
  				else
  					fig1 = figure('Color','white');
  					ax = axes('Parent',fig1);
  					res(k)=plot(ax,xx,yy,opts{:},'DisplayName',obj(k).label);
  				end
  			end
  		end
  		res(l==0)=[];%delete null handles
  	end
  end

end