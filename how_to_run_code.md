# Using CARC:
1. Navigate to your directory on your laptop so that when you type in `ls` you can see the right folder.
2. Use Cyberduck or `rsync` to upload your directory onto CARC by typing `rsync -vhatP myFolderName myUsername@wheeler.alliance.unm.edu:~/`
3. Type in your password.
4. Connect to the machine typing `ssh myUsername@wheeler.alliance.unm.edu`
5. Type in your password.
6. Navigate to where your shell script is on the machine.
7. Run your shell script by typing  `sh myScriptName.sh`. Example shell scripts are in the repository.
8. Once your code is done running, transfer a copy of your folder with results by typing `rsync -vhatP myUsername@wheeler.alliance.unm.edu:~/myFolderName newFolderName`

# Using `sftp` instead of `rsync` with Malarina as an example:
1. Navigate to your directory so that when you type in `ls` you can see the right folder.
2. Use `sftp` to upload your folder onto Malarina by typing `sftp remoteguest@malarina.unm.edu`
3. Type in the password.
4. Use the `put` command by typing `put -r myFolderName` and exit once it's transferred by typing `exit`
5. Connect to the machine using `ssh remoteguest@malarina.unm.edu`
5. Type in the password.
6. Navigate to where your shell script is on the machine.
7. Run your shell script by typing `sh myScriptName.sh`. Example shell scripts are in the repository.
- How to transfer a copy of your folder with results to your laptop:
8. Once your code is done running, navigate to the directory where you want your results to be.
9. Repeat steps 2-5, except replace "put" with "get".
