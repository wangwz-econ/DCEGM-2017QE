classdef cls_line
% Each object in the class cls_line is essentially a function or a line denoted by x and y values.
% This class defines certain methods to deal with a line.
% The most important methods are cls_upper_envlp() and cls_outer_refinement(). 
%   These two methods implement the grid refinement process in the paper 'The Endogenous Grid Method for
%   Discrete-Continuous Dynamic Choice Models with (or without) Taste Shocks (QE, 2017)'.

% Three properties:
%   cls_line.x, N by 1 vector
%   cls_line.y, N by 1 vector
%   label, character
% Together, cls_line.x and cls_line.y form a function (or a line).

% Methods: cls_len() to return the number of points
%					 cls_sort() to sort the points
%					 cls_interp() to compute the function
%					 cls_insert() to insert a new point
%					 cls_chop() to break one line into two
% 				 cls_grow() to join multiple lines together
%					 cls_thinout() to delete points by index
%					 cls_diff() to compare two polylines
%					 cls_upper_envlp() to find upper envelope of a set of lines
%					 cls_outer_refinement() to refine the endogenous wealth grid

properties (Access=public)
  x = [];    % x-grid points
  y = [];    % function values
  label = "" % label of this line
end

methods

  % cls_line %----------------------------------------------
  function obj = cls_line(varargin)
    % Initialize the object.
    % ! x and y must be vectors of the same size. After initialization, obj.x and obj.y are both column vectors.
    if nargin > 0 
      obj.x = reshape(varargin{1}, [], 1); % reshape x as a column vector
    end
    if nargin > 1
      obj.y = reshape(varargin{2}, [], 1); % reshape y as a column vector
    end
    if nargin > 2 && ischar(varargin{3})
      obj.label = varargin{3};
    end
  end

  % cls_len %----------------------------------------------
  function res = cls_len(obj)
    % Number of grid points in this line/function.
    % Also can check the equality of length of obj.x and obj.y.
    % ! obj can be a scalar or array. 
    % ! When obj is an array, res is an array of the same size with entry k storing the length of x and y of obj(k).
    res = NaN(numel(obj), 1);
    for k = 1:numel(obj)
      res(k) = numel(obj(k).x);
      if res(k)~=numel(obj(k).y)
        error('The %d th obj is not a line (mismatch in number of points): x has %d while y has %d elements', ...
          k, numel(obj(k).x), numel(obj(k).y))
      end
    end % end the loop for all elements in the obj array
    res = reshape(res, size(obj));
  end

  % cls_sort %----------------------------------------------
  function res = cls_sort(obj)
    % Sort the obj.x vector to better represent a function.
    % Return the sorted line (another cls_line object or object array), original obj is not modified.
    % ! obj can be a scalar or array.
    for i = 1:numel(obj)
      [res_x, sort_index] = sort(obj(i).x);
      res_y = obj(i).y(sort_index);
      res(i) = cls_line(res_x, res_y, obj(i).label);
    end % end the loop for all elements in the obj array
    res = reshape(res, size(obj));
  end % end the cls_sort function

  % cls_interp %----------------------------------------------
  function [res, extrapflag] = cls_interp(obj, vx)
    % Compute the linearly interploated/extrapolated values in given points .
    % ! vx is a row or column vector.
    % ! If obj is a cls_line scalar:
    % !   res is a row vector of the same size of vx, storing the interpolated/extrapolated values of vx.
    % !   extrapflag is a logical row vector of the same size of vx (1 means this value is extrapolated).
    % ! If obj is a cls_line array:
    % !   res is a cell of the same size as the obj array, with entry k containing the
    % !     interploated/extrapolated values (as a row vector) for the corresponding function obj(k).
    % TODO obj(k).x must only unique points but need not to be sorted.
    vx = reshape(vx, 1, []);
    if numel(obj)==1
      try
        res = interp1(obj.x, obj.y, vx, "linear", "extrap");
      catch ME
        steps = reshape(obj.x(2:end)-obj.x(1:end-1),[],1);
        fprintf('\nInterpolation error!\n\nRepeated points: %d\nInverted steps: %d\n', ...
               sum(abs(steps)<eps), sum(steps<0));
        fprintf('Identifier: %s\nCall stack: ', ME.identifier);
      end % end the try-catch statement
      extrapflag = ~(vx>=min(obj.x) & vx<=max(obj.x)); % returns true if extrap is needed
    else
      res = cell(size(obj));
      extrapflag = cell(size(obj));
      for k=1:numel(obj)
        try
          res{k} = interp1(obj(k).x, obj(k).y, vx, 'linear', 'extrap');
        catch ME
          steps = reshape(obj(k).x(2:end)-obj(k).x(1:end-1), [], 1);
          fprintf('\n\nInterpolation error!\nLine Number=%d,\nRepeated points: %d\nInverted steps: %d\n', ...
                 k, sum(abs(steps(k))<eps), sum(steps(k)<0) );
          fprintf('Identifier: %s\nCall stack: ',ME.identifier);
        end % end the try-catch statement
        extrapflag{k} = ~( vx>=min(obj(k).x) & vx<=max(obj(k).x) );
      end % end the for loop to go through all entries in the object array
    end % end the if-else statement to check the number of elements in the obj input
  end % end the cls_interp function

  % cls_grow %----------------------------------------------
  function res = cls_grow(obj, in2, front)
    % Add another line segement in2 to the end (or the front) of ojb.
    % ! If obj is a cls_line scalar:
    % !   in2 is also a cls_line scalar.
    % !   If front = false (default), then add in2 to the end of obj; otherwise, to the front of obj.
    % !   Return a cls_line object scalar res.
    % ! If obj is a cls_line array:
    % !   in2 is also a cls_line array of the same size of obj.
    % !   If front = false (default), then add in2 to the end of obj; otherwise, to the front of obj. 
    % !   Return a cls_line object array res. (use res(index) to access certain lines).
    if nargin<3
      front = false;
    end
    for k1 = 1:numel(obj)
      if front
        res(k1) = cls_line( ...
          [reshape(in2(k1).x, 1, []), reshape(obj(k1).x, 1, [])], ...
          [reshape(in2(k1).y, 1, []), reshape(obj(k1).y, 1, [])], ...
          sprintf('%s + %s', in2(k1).label, obj(k1).label));
      else
        res(k1) = cls_line( ...
          [reshape(obj(k1).x,1,[]), reshape(in2(k1).x, 1, [])], ...
          [reshape(obj(k1).y,1,[]), reshape(in2(k1).y, 1, [])], ...
          sprintf('%s + %s',obj(k1).label,in2(k1).label));
      end % end the if statement to check front argument
    end % end the loop to go through all elements in the obj array
    res = reshape(res, size(obj));
  end % end the cls_grow function

  % cls_insert %----------------------------------------------
  function res = cls_insert(obj, x, y, j)
    % Insert points (x, y) after position j in the current line obj.
    % ! j=0 means to add the points to the front, j>=obj.len means to add the points to the end.
    % ! If j is missing, the point is added to the end.
    % ! x, y can only be vectors (row or column vectors of the same size) or scalars.
    % ! obj can be a cls_line scalar or array. If it is an array, then insert points (x,y) into each cls_line object.
    res = obj;
    if numel(x)>1 % here x and y are two vectors
      mx = numel(x);
      my = numel(y);
      if mx ~= my
        warning('The x and y inputs do not have the same length!')
      end
      x = reshape(x,[],1); y = reshape(y,[],1);
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
    elseif numel(obj)>1 % here the (x,y) are two scalars, obj is a cls_line object array
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
    else % here the (x,y) are two scalars, obj is a cls_line object scalar
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
  end % end the cls_insert function

  % cls_thinout %----------------------------------------------
  function res = cls_thinout(obj, index)
    % Remove the indexed points from obj.
    % ! index should be a scalar or a vector (row or column vector).
    % ! obj can be a cls_line object scalar or array. If it is an array, removing the indexed points in all entries.
    res = obj;
    for k = 1:numel(res)
      ii = intersect(1:numel(res(k).x), index);
      remove_index = reshape(ii, 1, []);
      res(k).x(remove_index) = [];
      res(k).y(remove_index) = [];
    end
  end % end the cls_thinout function

  % cls_diff %----------------------------------------------
  function indx = cls_diff(obj, pl2, significance)
    % Get the indexes of points in obj that are not in pl2.
    % Note: differences in either x or y will both be stored.
    % ! pl2 should be a cls_line object scalar.
    % ! obj can be a cls_line object scalar or array.
    % ! If obj is a cls_line object array, indx is a cell; if obj is a cls_line object scalar, indx is a vector.
    if nargin < 3
      significance = 5; % equality is measured up to 10^-signif
    end
    x1 = round(pl2(1).x * (10^significance)) * 10^(-significance);
    y1 = round(pl2(1).y * (10^significance)) * 10^(-significance);
    for k = 1:numel(obj)
      a = round(obj(k).x * (10^significance)) * 10^(-significance);
      b = round(obj(k).y * (10^significance)) * 10^(-significance);
      indx{k} = find(~ismember(a, x1) | ~ismember(b, y1));
    end
    if numel(indx)==1
      indx = indx{1}; 
    else
      indx = reshape(indx, size(obj));
    end
  end % end the cls_diff function

  % cls_chop %----------------------------------------------
  function [res1, res2] = cls_chop(obj, j, repeat)
    % This function separates the grid into two parts at the index j position.
    % ! j must be a scalar.
    % ! obj can be a cls_line object scalar or array. 
    % ! res1 and res2 are cls_line object arrays of the same size as obj.
    % If repeat=true, the boundary points are repeated in both resulting cls_line objects, i.e., 
    %   the two parts are 1,..., j and j, ..., end.
    % If repeat=false, the two parts are 1,..., j and j+1, ..., end.
    if nargin<2
      error('Have to have one input for .cls_chop')
    end
    for k=1:numel(obj) % chop each polyline if there are many
      if j>obj(k).cls_len
        warning('Producing empty polyline by chopping at index j>len')
        j = obj(k).len;
      end
      res1(k) = cls_line(obj(k).x(1:j), obj(k).y(1:j), sprintf('%s (1:%d)', obj(k).label, j));
      if nargout > 1
        if repeat
          res2(k) = cls_line(obj(k).x(j:end), obj(k).y(j:end), ...
            sprintf('%s (%d:%d)', obj(k).label, j, obj(k).cls_len));
        else
          res2(k) = cls_line(obj(k).x(j+1:end), obj(k).y(j+1:end), ...
            sprintf('%s (%d:%d)', obj(k).label, j+1, obj(k).cls_len));
        end % end the if-else statement to check whether repeat==true
      end % end the if-else statement to check the number of output arguments
    end % end the for loop to go through all elements in obj array
    res1 = reshape(res1, size(obj));
    if nargout > 1
      res2 = reshape(res2, size(obj));
    end
  end % end the cls_chop function

  % cls_upper_envlp %----------------------------------------------
  function [res, intersections] = cls_upper_envlp(obj, fullinterval)
    % This function computes the upper envelope over the array of cls_line objects.
    % ! obj must be a cls_line object column vector of the same length in the sense of obj(k).cls_len!
    % TODO All entries in obj should be sorted (obj(k).x are unique and increasing), and are treated as sorted.
    % ! If fullinterval=false, res is a cls_line object containing the upper envelope of only the overlapping segments.
    % ! If fullinterval=true, res is a cls_line object containing the upper envelope of the whole interval (union of obj(k).x).

    if numel(obj)==1
      warning('Upper envelope is meant for an array of polylines')
      res = obj;
      if nargout>1
        intersections = cls_line([], [], "intersection points");
      end
      return
    end % end the if statement to check the number of elements in obj array

    % step 1: identify the query points.
    l = obj.cls_len;            % check that x and y are of the same size
    obj = obj(l>0);             % disregard all polylines of zero length
    pt = obj(1).x;              % collect all the x points in sorted row vector
    for union_k = 2:numel(obj)
      pt = union(pt, obj(union_k).x);
    end
    pt = sort(pt);  
    % ! When obj is a cls_polyline array, ojb.x returns all obj(k).x as vectors of possible different lengths.

    % step 2: interpolate all lines on all points recording the extrapolation cases
    [intr, extrflag] = obj.cls_interp(pt); % intr and extrflag are size(obj) by 1 cells, each entry containing a row vector of interpolated values
    intr = cell2mat(intr); extrflag = cell2mat(extrflag); % now they become size(obj) by size(pt) matrices
    if ~fullinterval
      % disregard points where at least one line is extrapolated
      mask = (sum(extrflag, 1) > 0); % sum over rows (different cls_line objects), mask is a logical row vector
      intr(:, mask) = [];
      pt(mask) = [];
      n = sum(~mask); % number of points in the overlap region
    else
      % disregard only points where particular lines are extrapolated (full interval!)
      intr(extrflag) = -Inf; % -Inf so that the corresponding extrapolated values won't construct the upper envelope
      n = numel(pt);
    end

    % step 3: find lines on the top
    maxintr_vec = max(intr); % maximization across rows (across different cls_line objects), exactly what we want
    maxintr = repmat(maxintr_vec, size(intr,1), 1); % a matrix of the same size as intr, rows are the same (as the upper envelope)
    top = (intr == maxintr);

    % step 4.1: build up the upper envelope for the first point
    res = cls_line(pt(1), maxintr(1,1), 'upper envelope');
    if nargout>1
      intersections = cls_line([], [], "intersection points");
    end
    k0 = find(top(:,1), 1, 'first'); % index of top line

    % step 4.2: go through all points
    for j = 2:n
      k1 = find(top(:,j), 1, "first"); % index of next top line
      if k1 ~= k0 % this case means that the upper envelope changes line segments for k0 to k1
        
        % check if there is an intersection point between the lines:
        ln1=k0; ln2=k1;             % check if there is an intersection between lines ln1 and ln2
        pt1 = pt(j-1); pt2 = pt(j); % possible intersection lies in interval [pt1, pt2]
        [y1, extr1] = obj(ln1).cls_interp([pt1 pt2]);
        [y2, extr2] = obj(ln2).cls_interp([pt1 pt2]); % and these function values (maybe extrapolated)

        % check that neither is extrapolated in both points
        % and that intersection point is inside the inverval <= function values are different at the borders
        if all(~[extr1 extr2]) && all(abs(y1-y2)>0)
          % find the intersection point or points
          while true
            fun_forintersection = @(x) obj(ln2).cls_interp(x) - obj(ln1).cls_interp(x);
            pt3 = fzero( fun_forintersection, [pt1 pt2] );
            pt3f = obj(ln1).cls_interp(pt3);
            % check if there are lines above the found intersection
            [intr2, exrt2] = obj.cls_interp(pt3); % interpolate all lines in the new point, return size(obj) by 1 cells
            intr2 = cell2mat(intr2); exrt2 = cell2mat(exrt2); % now they become size(obj) by 1 vectors
            intr2(exrt2) = -Inf; % disregard the extrapolated points
            maxintr2 = repmat(max(intr2), size(intr2, 1), 1);
            ln3 = find(intr2==maxintr2, 1, 'first');
            if ln3==ln1 || ln3==ln2 % there are no other functions above!
              % add the intersection point
              res = res.cls_insert(pt3, pt3f, res.cls_len); % add the intersection point to the end
              if nargout > 1
                intersections = intersections.cls_insert(pt3, pt3f);
              end
              if ln2 == k1 % maybe there are some more intersections before next point?
                % no actually, because the left line is the one we started with
                break
              else % indeed, so update the interval of new search
                ln1 = ln2;
                pt1 = pt3;
                ln2 = k1;
                pt2 = pt(j);
              end
            else % there is another function above ln1 and ln2 at calculated intersection point pt3
              % so, it is not on the upper envelope. need to search again.
              ln2 = ln3; % new candidate
              pt2 = pt3;
            end % end the if statement to check whether the upper envelope ln3 is either ln2 or ln1
          end % end the while loop to get all the intersections
        end % end the if statement to check there is an intersection when changing line segments of the upper envelope
      end % end the if statement to check line segments of the upper envelope have changed
      if ismember(pt(j), obj(k1).x) || j==n
        res = res.cls_insert(pt(j), maxintr(1, j)); % add new grid point to the end
        res.x = reshape(res.x, [], 1); % since we add points one by one, 
        res.y = reshape(res.y, [], 1); %   this step ensures that all points are stacked in a column vector
      end
      k0 = k1; % next step
    end % end the for loop to go through all points in the pt vector
  end % end the cls_upper_envlp function

  % cls_outer_refinement %----------------------------------------------
  function [res, indxremoved, newdots] = cls_outer_refinement(obj)
		% Conduct the endogenous wealth grid refinement as in the paper 'The Endogenous Grid Method for Discrete-Continuous 
    %   Dynamic Choice Models with (or without) Taste Shocks, Quantitative Economics, 2017'.
    % TODO The first step is to detect the non-monotonic (non-increasing) region in obj.x, 
    % TODO   and then chop the obj into several pieces, with each piece the x values are monotonic.
    % TODO The second step is to construct the upper envelope of the lines in different segments using cls_upper_envlp.
    % TODO The inner process of constructinng upper envelope is actually refinement to obj.x.
    % ! obj can be cls_line scalar or array. This function is imposed elementwise.
		% ! res is a cls_line scalar or array after refining obj.
    % ! indxremoved is a cell or a vector of the same size as obj. Entry k stores indexes of deleted points in obj(k).
    % ! newdots is a cls_line scalar or array of size as size(obj). Entry k stores new points from line intersections.

		for k=1:numel(obj)
			cur = obj(k);
			% identify the loop-back regions
			ii = (cur.x(2:end)>cur.x(1:end-1)); % zeros for the loopback regions (due to non-increasing in cur.x)
			sect = cls_line; 
			i = 1;
			while true
				j = find(ii~=ii(1), 1, 'first');
				if isempty(j)	% exit the loop if there are no more loop-backs (non-increasing regions in cur.x)
					if i>1	% if sections already started, add the last one
						sect(i) = cur;
					end
					break;
				end
				[sect(i), cur] = cls_chop(cur, j, true); % true to repeat the boundary points
				ii(1:j-1) = []; % chop ii array accordingly
				i = i + 1;
			end % end the while loop to go through all loopback regions
			% perform secondary envelope if sections created
			if numel(sect)>1
				sect = sect.cls_sort; % sort all sections since half of them are in opposite direction
        sect = reshape(sect, [], 1);
				[res(k), newdots(k)] = sect.cls_upper_envlp(true); %true for full interval
				% removed points indexes
				indxremoved{k} = cls_diff(obj(k), res(k), 10); % index of dots in obj(k) but not in res(k)
			else
				% without loopbacks -- just copy the input
				res(k)=obj(k);
				indxremoved{k}=[];
				newdots(k) = cls_line;
			end
		end % end the for loop to go through all obj entries
		% return in same dimensions
		res = reshape(res, size(obj));
		indxremoved = reshape(indxremoved, size(obj));
		if numel(obj)==1
			indxremoved=indxremoved{1}; % for a single polyline return indx as a vector
		end
	end % end the cls_outer_refinement function

end % end method

end % end classdef
