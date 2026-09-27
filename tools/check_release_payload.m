function report=check_release_payload(outputDir,root)
%CHECK_RELEASE_PAYLOAD Verify ZIP/MLTBX contents without installing an add-on.
%   Compares the full released payload byte-for-byte with ROOT, verifies
%   the stable identifier, version/name and exactly two installation paths.
    if nargin<2, root=fileparts(fileparts(mfilename('fullpath'))); end
    scratch=tempname; mkdir(scratch); cleanup=onCleanup(@() rmdir(scratch,'s'));
    zipRoot=fullfile(scratch,'zip'); installerRoot=fullfile(scratch,'installer');
    unzip(fullfile(outputDir,'ContourRoots.zip'),zipRoot);
    unzip(fullfile(outputDir,'ContourRoots.mltbx'),installerRoot);
    payload=fullfile(installerRoot,'fsroot');
    expected={};
    for item={'matlab','examples','docs','setup_contourroots.m','README.md','LICENSE', ...
            'THIRD_PARTY_NOTICES.md','CITATION.cff','CHANGELOG.md','CONTRIBUTING.md','ROADMAP.md'}
        if isfolder(fullfile(root,item{1}))
            d=dir(fullfile(root,item{1},'**','*')); d=d(~[d.isdir]);
            entries=cellfun(@(p) strrep(p,[root filesep],''),fullfile({d.folder},{d.name}),'UniformOutput',false);
        else
            entries=item;
        end
        expected=[expected entries]; %#ok<AGROW>
    end
    expected=strrep(expected,'\','/');
    expected=sort(expected(~startsWith(expected,'docs/development/') & ...
        ~endsWith(expected,'.DS_Store') & ~endsWith(expected,'.asv')));
    assert(isequal(expected,files(zipRoot)),'check_release_payload:ZipFiles','ZIP payload differs from release manifest.');
    assert(isequal(expected,files(payload)),'check_release_payload:InstallerFiles','MLTBX payload differs from ZIP.');
    for k=1:numel(expected)
        original=bytes(fullfile(root,expected{k}));
        assert(isequal(original,bytes(fullfile(zipRoot,expected{k}))), ...
            'check_release_payload:ZipBytes','ZIP changed %s.',expected{k});
        assert(isequal(original,bytes(fullfile(payload,expected{k}))), ...
            'check_release_payload:InstallerBytes','MLTBX changed %s.',expected{k});
    end
    configuration=xmlread(fullfile(installerRoot,'metadata','configuration.xml'));
    nodes=configuration.getElementsByTagName('matlabPath'); paths=cell(1,nodes.getLength);
    for k=1:numel(paths), paths{k}=char(nodes.item(k-1).getTextContent); end
    assert(isequal(sort(paths),{'/matlab','/matlab/delay'}), ...
        'check_release_payload:Path','Installer path metadata changed.');
    addon=fileread(fullfile(installerRoot,'metadata','addonProperties.xml'));
    assert(contains(addon,'6f1b2e7c-3c1a-4a55-9b7e-0c4b8f2d9a31'),'check_release_payload:ID','Installer identity changed.');
    core=fileread(fullfile(installerRoot,'metadata','coreProperties.xml'));
    source=fileread(fullfile(root,'matlab','contourroots_version.m'));
    v=char(regexp(source,'''(\d+\.\d+\.\d+)''','tokens','once'));
    assert(contains(core,['<cp:version>' v '</cp:version>'])&&contains(core,'<dc:title>ContourRoots</dc:title>'), ...
        'check_release_payload:Version','Installer version/name differs from source.');
    report=struct('files',numel(expected),'byteIdentical',true,'installationPaths',{paths},'version',v);
    fprintf('Verified %d byte-identical files in ZIP and MLTBX; installer metadata retained.\n',numel(expected));
end
function names=files(root)
    d=dir(fullfile(root,'**','*')); d=d(~[d.isdir]);
    names=sort(strrep(strrep(fullfile({d.folder},{d.name}),[root filesep],''),'\','/'));
end
function b=bytes(file)
    fid=fopen(file,'rb'); assert(fid>=0,'check_release_payload:Read','Cannot read %s.',file);
    cleanup=onCleanup(@() fclose(fid)); b=fread(fid,Inf,'*uint8');
end
