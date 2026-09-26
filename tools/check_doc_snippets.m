function report = check_doc_snippets(files)
%CHECK_DOC_SNIPPETS Execute the MATLAB code blocks of Markdown documents.
%   CHECK_DOC_SNIPPETS(FILES) extracts every fenced code block marked
%   ```matlab from each Markdown file in FILES (char, string or cellstr)
%   and runs the blocks of one file in order, in a shared workspace, like
%   the sections of a live script. Figures are created invisibly.
%
%   A block is skipped when the line immediately before its opening fence
%   is the HTML comment <!-- no-test --> (for example, download commands).
%
%   REPORT = CHECK_DOC_SNIPPETS(...) returns a table with the file, block
%   number, first line of the block, and status. An error is thrown after
%   all files are processed if any block failed.

    files = cellstr(files);
    rows = cell(0,4);
    failures = 0;
    oldVisible = get(groot,'defaultFigureVisible');
    before = findall(groot,'Type','figure');
    set(groot,'defaultFigureVisible','off');
    restore = onCleanup(@() cleanup(oldVisible,before));
    for i = 1:numel(files)
        blocks = extract_blocks(files{i});
        [status,messages] = run_blocks(blocks);
        for b = 1:numel(blocks)
            first = strtrim(strtok(blocks(b).code,newline));
            rows(end+1,:) = {string(files{i}), b, string(first), string(status{b})}; %#ok<AGROW>
            if strcmp(status{b},'failed')
                failures = failures + 1;
                fprintf(2,'FAILED %s, block %d (line %d):\n%s\n\n', ...
                    files{i}, b, blocks(b).line, messages{b});
            end
        end
    end
    report = cell2table(rows,'VariableNames',{'File','Block','FirstLine','Status'});
    if failures > 0
        error('check_doc_snippets:Failed','%d documentation block(s) failed.',failures);
    end
end

function blocks = extract_blocks(file)
    lines = splitlines(string(fileread(file)));
    blocks = struct('code',{},'line',{},'skip',{});
    k = 1;
    while k <= numel(lines)
        if strtrim(lines(k)) == "```matlab"
            skip = k > 1 && contains(lines(k-1),'<!-- no-test -->');
            start = k + 1;
            k = k + 1;
            while k <= numel(lines) && strtrim(lines(k)) ~= "```"
                k = k + 1;
            end
            code = strjoin(lines(start:k-1), newline);
            blocks(end+1) = struct('code',char(code),'line',start-1,'skip',skip); %#ok<AGROW>
        end
        k = k + 1;
    end
end

function [cr__status,cr__messages] = run_blocks(cr__blocks)
    % All blocks of one document run in THIS function's workspace, so that
    % variables defined in one block are visible in the next. Local names
    % use a cr__ prefix to avoid clashing with names used in the snippets.
    cr__status = repmat({'passed'},1,numel(cr__blocks));
    cr__messages = repmat({''},1,numel(cr__blocks));
    cr__folder = tempname; mkdir(cr__folder);
    cr__remove = onCleanup(@() rmdir(cr__folder,'s'));
    for cr__b = 1:numel(cr__blocks)
        if cr__blocks(cr__b).skip
            cr__status{cr__b} = 'skipped'; continue
        end
        cr__file = fullfile(cr__folder,sprintf('cr_doc_block_%d.m',cr__b));
        cr__fid = fopen(cr__file,'w');
        fprintf(cr__fid,'%s\n',cr__blocks(cr__b).code); fclose(cr__fid);
        try
            evalc('run(cr__file)');
        catch cr__err
            cr__status{cr__b} = 'failed';
            cr__messages{cr__b} = getReport(cr__err,'basic');
        end
    end
end

function cleanup(oldVisible,before)
    % Close only the figures created by the snippets.
    close(setdiff(findall(groot,'Type','figure'),before));
    set(groot,'defaultFigureVisible',oldVisible);
end
