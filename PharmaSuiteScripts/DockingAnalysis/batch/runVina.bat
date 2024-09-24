@echo off
for %%a in ("%~dp0\*.pdbqt) do (
"C:\adt\Lib\site-packages\AutoDockTools\vina.exe" --ligand "%%a" --config configNew.txt
) 

pause