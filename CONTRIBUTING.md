# Contributing to the database: 
We anticipate that ForC-db will be of great value for various initiatives seeking to understand and manage the role of tropical forests in the global C cycle, and we encourage authors to contribute data. We hope that this will provide an avenue for data contributors to collaborate on projects using the data that they provide, and we encourage database users to inform researchers who contribute to the database of intended use of the data and to discuss potential collaboration (see [`DATA USE POLICY.md`](https://github.com/forc-db/ForC/blob/master/DATA%20USE%20POLICY.md)). However, the database is open access, and data contributors assume the risk that some users may not adhere to our best practices guidelines. 

Instructions for formally adding data to the master branch of ForC-db are given below.

Alternatively, if you want to post data that has not been formatted for ForC-db, you can do so by creating a [Github issue](https://help.github.com/articles/creating-an-issue). That way the existence of a dataset can be recorded and discussed without requiring anyone to format it until someone is motivated to do so (e.g., because it will help them with their research).


## Basic instructions for editing ForC-db: 
1.	Download latest version of files.
2.	Check if data has already been entered in database:  
    a.	In **measurements**:  
      i.  Check *citations.doi*, *citations.author*, *citations.year*, *citations.title*  
      ii. Check *sites.sitename*  
      iii.Check *plot.name*  
    b.	In **sites**:  
      i.	Check *sites.sitename*, *lat*, *lon*, *masl*   
    c.	If the site already exists in the database and site-related data matches, please use the same *sites.sitename*  
    d.	If the **plot.name** already exists, please use the same plot.name  
3.	If you are planning to make a substantial contribution that will take more than a few days to complete, we recommend creation of an **issue** in Github to alert database owners about the changes intended to be made. This will avoid potential duplication of efforts and allow database owners to provide any necessary guidance.  
4.	Create a new **branch**, including your username in the branch name.  
5.	Modify files, ensuring that  
    a.	If entering new data, begin entering data in the primary tables. We recommend to start with **measurements**, followed by **sites**, **plothistory**, **methodology** and **allometry** (if applicable).  
    b.	The remaining tables (**variables**, **disttype**, **pft**) will require updates only when new variables, disturbance type, plant functional type are introduced. These should rarely be updated.  
    c.	Refer to [readme_ForC.pdf](ForC/readme_ForC.pdf) for details on fields in each table.   
    d.	Every cell should be populated. For missing data, please use missing data codes (Additional Table 1 in [readme_ForC.pdf](ForC/readme_ForC.pdf)) to indicate reason for missing values.  
    e.	Conduct basic QA/QC to make sure additions/changes are made in all linked files.  
6.	Add new rows to the bottom of existing files (or replace old files with updated files), and commit changes to the newly created branch.  
7.	Open a [pull request](https://help.github.com/articles/creating-a-pull-request/). This will compare changes across the master branch and the new branch.  
8.	Database authors will approve pull request and merge changes. Include the words “fixes #__ <insert Issue number> to automatically close issue. Otherwise, manually close issue and refer to pull request. Delete branch when changes have been merged.  
