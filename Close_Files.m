fclose(FNOUT);
if (EO_PRINT)
    fclose(FNEVEN); fclose(FNODD);
end
if (THERMO)
    fclose(FNTHERM); fclose(FNTHERMAVG);
end
if (OPMAT)
    fclose(FNOPMAT);
end