

set "python_path=C:\adt\python.exe"
set "script_path=C:\adt\Lib\site-packages\AutoDockTools\prepare_ligand4.py"

for %%a in ("%~dp0\LigandEmmet*.pdb") do (
    %python_path% %script_path% -l "%%a" -A bonds -v -U nphu 
)
pause