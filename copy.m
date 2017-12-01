% Copy function - replacement for matlab.mixin.Copyable.copy() to create object copies
function newObj = copy(obj)
    try
        % R2010b or newer - directly in memory (faster)
        objByteArray = getByteStreamFromArray(obj);
        newObj = getArrayFromByteStream(objByteArray);
    catch
        % R2010a or earlier - serialize via temp file (slower)
        fname = [tempname '.mat'];
        save(fname, 'obj');
        newObj = load(fname);
        newObj = newObj.obj;
        delete(fname);
    end
end