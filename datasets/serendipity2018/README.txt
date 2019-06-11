Summary
=======

We conducted an experiment in MovieLens (http://movielens.org), where we asked users how serendipitous particular movies were to them. This dataset contains user answers to our questions and additional information, such as past ratings of these users, recommendations they received before replying to our survey and movie descriptions. The dataset is meant for research regarding serendipity in recommender systems, such as analysis of serendipitous movies or offline evaluation of serendipity-oriented recommendation algorithms.
The dataset was generated on January 15, 2018. The data are contained in the files ` answers.csv`, ` movies.csv`, ` recommendations.csv`, ` tag_genome.csv`, ` tags.csv` and ` training.csv`. Overall, there are 10,000,000 ratings (2,150 ratings stored in `answers.csv` and 9,997,850 in `training.csv`). More details about the contents and use of all these files follow.
Note that this dataset is slightly different from the one used in the linked publication (section `Citation`), as these data were updated, while we were reporting the results.
In this dataset, each selected user had rated at least 1 movie. No demographic information is included. Each user is represented by an id, and no other information is provided.

This and other GroupLens datasets are publicly available for download at <http://grouplens.org/datasets/>.

Usage License
=============

Neither the University of Minnesota nor the University of Jyväskylä nor any of the researchers involved can guarantee the correctness of the data, its suitability for any particular purpose, or the validity of results based on the use of the dataset. The dataset may be used for any research purposes under the following conditions:

* The user may not state or imply any endorsement from the University of Jyväskylä, the University of Minnesota or the GroupLens Research Group.
* The user must acknowledge the use of the dataset in publications resulting from the use of the dataset (see below for citation information).
* The user may not redistribute the data without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from a faculty member of the GroupLens Research Project at the University of Minnesota.
* The executable software scripts are provided "as is" without warranty of any kind, either expressed or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. The entire risk as to the quality and performance of them is with you. Should the program prove defective, you assume the cost of all necessary servicing, repair or correction.

In no event shall the University of Minnesota, the University of Jyväskylä and their affiliates or employees be liable to you for any damages arising out of the use or inability to use these programs (including but not limited to loss of data or data being rendered inaccurate).

If you have any further questions or comments, please email <grouplens-info@umn.edu>


Citation
========

To acknowledge use of the dataset in publications, please cite the following paper:
Denis Kotkov, Joseph A. Konstan, Qian Zhao, and Jari Veijalainen. 2018. Investigating Serendipity in Recommender Systems Based on Real User Feedback. In Proceedings of SAC 2018: Symposium on Applied Computing , Pau, France, April 9–13, 2018 (SAC 2018), 10 pages. DOI: 10.1145/3167132.3167276

Further Information About GroupLens
===================================

GroupLens is a research group in the Department of Computer Science and Engineering at the University of Minnesota. Since its inception in 1992, GroupLens's research projects have explored a variety of fields including:

* recommender systems
* online communities
* mobile and ubiquitious technologies
* digital libraries
* local geographic information systems

GroupLens Research operates a movie recommender based on collaborative filtering, MovieLens, which is the source of these data. We encourage you to visit <http://movielens.org> to try it out! If you have exciting ideas for experimental work to conduct on MovieLens, send us an email at <grouplens-info@cs.umn.edu> - we are always interested in working with external collaborators.


Content and Use of Files
========================

Formatting and Encoding
-----------------------

The dataset files are written as [comma-separated values](http://en.wikipedia.org/wiki/Comma-separated_values) files with a single header row. Since some values contain commas (`,`) in movies.csv, all the values are enclosed in double-quotes (`"`) in this file. These files are encoded as UTF-8. If accented characters in movie titles or tag values (e.g. Misérables, Les (1995)) display incorrectly, make sure that any program reading the data, such as a text editor, terminal or script, is configured for UTF-8.

User Ids
--------

Their ids have been anonymized. User ids are consistent between `answers.csv`, `recommendations.csv`, `tags.csv` and `training.csv` (i.e. the same id refers to the same user across the four files).


Movie Ids
---------

Movie ids are consistent with those in other MovieLens datasets and those used on the MovieLens web site (e.g., id `1` corresponds to the URL <https://movielens.org/movies/1>). Movie ids are consistent between ` answers.csv`, ` recommendations.csv`, `movies.csv`, `tag_genome.csv`, `tags.csv` and ` training.csv` (i.e. the same id refers to the same movie across these six data files).

Timestamps
---------

A timestamp indicates date and time when a user rated a particular movie. Timestamps represent seconds since midnight Coordinated Universal Time (UTC) of January 1, 1970.

Ratings
---------

Users assign ratings to movies on the scale from 0.5 star till 5 stars with the granularity of 0.5 star. Ratings represent how much users enjoyed movies. The higher the rating the higher the enjoyment.

Answer Data File Structure (answers.csv)
-----------------------------------------

To investigate serendipity, we invited users to complete our survey on April 1, 2017. We selected users who rated at least five movies with a rating of at least 3.5 stars from December 30, 2016 till March 30, 2017. Among users who joined MovieLens after November 30, we selected those, who rated at least five movies with a rating of at least 3.5 stars one month after their registration. In our survey, we asked users to answer a question and rate forty statements about each of five movies we picked for these users (five questions and forty statements, overall). For each user, we picked five the least popular movies (movies with the fewest numbers of ratings in MovieLens) among movies this user gave at least 3.5 stars during the three months or less for users who recently joined the system. Not all the users rated all the statements and answered all the questions about each movie.

We invited users to our survey on April 1, 2017 and generated the dataset for the linked publication in May, 2017. However, a few more users replied to our survey and some users updated their survey replies and ratings after we generated the dataset. This dataset contains the updated data. There are 2,150 records. The file contains the following fields:

userId – user id (481 users)
movieId – movie id (1,678 movies)
rating – see the description above
timestamp – timestamp, which indicates when the user gave the rating. The ratings in this file were given or updated from Decemer 30, 2016 till January 6, 2018
predictedRating – the rating, which MovieLens predicted for the user and the movie just before the user rated the movie. In this file, 20 ratings have missing predicted rating
The following fields correspond to user ratings of the statements in our survey. Ratings are given using the scale, where 1 corresponds to `strongly disagree`, 2 – `disagree`, 3 – `neither agree nor disagree`, 4 – `agree`, 5 – `strongly agree`, NA – `don't remember`
s1 – `The first time I heard of this movie was when MovieLens suggested it to me.`
s2 – `MovieLens influenced my decision to watch this movie.`
s3 – `I expected to enjoy this movie before watching it for the first time.`
s4 – `This is the type of movie I would not normally discover on my own; I need a recommender system like MovieLens to find movies like this one.`
s5 – `This movie is different (e.g., in style, genre, topic) from the movies I usually watch.`
s6 – `I was (or, would have been) surprised that MovieLens picked this movie to recommend to me.``
s7 – `I am glad I watched this movie.`
s8 – `Watching this movie broadened my preferences. Now I am interested in a wider selection of movies.`
q – answer to the question `When did you watch this movie for the first time?`, where 1 corresponds to `past week`, 2 – `past month`, 3 – `1-6 months`, 4 – `6-12 months`, 5 – `1-3 years`, 6 – `>3 years ago`, 7 – `don't remember`.
The following binary fields are calculated based on ratings of the statements above (s1, s2, s3, s4, s5, s6, s7 and s8). For example, s_ser_find is TRUE if s1 > 3 and s4 > 3, and FALSE otherwise. For more details, please see the linked publication.
s_ser_rel – strict serendipity (relevant)
s_ser_find – strict serendipity (find)
s_ser_imp – strict serendipity (implicit)
s_ser_rec – strict serendipity (recommend)
m_ser_rel – motivational serendipity (relevant)
m_ser_find – motivational serendipity (find)
m_ser_imp – motivational serendipity (implicit)
m_ser_rec – motivational serendipity (recommend)

Rating Data File Structure (training.csv)
-----------------------------------------

This file contains all user ratings given or updated from November 11, 2009 till January 6, 2018 excluding ratings in answers.csv, but including at least one rating of each user from answers.csv. Overall, training.csv contains 9,997,850 ratings with the following fields:
userId – user id (104,661 users)
movieId – movie id (49,151 movies)
rating – see the description above
timestamp – timestamp, which indicates when the user gave the rating.

Recommendation Data File Structure (recommendations.csv)
-----------------------------------------

This file contains eight recommendations provided for each user in answers.csv before the survey during the period of time, from November 1, 2016 till May 1, 2017. This data is missing for the user with the id 249,371, who did not receive recommendations during the indicated period of time. The file contains the following fields:
userId – user id (480 users)
movieId – movie id (1,181 movies)

Movie Data File Structure (movies.csv)
-----------------------------------------

The file contains descriptions of the movies from `answers.csv`, `recommendations.csv` and `training.csv`. Since some values contain commas (`,`), all the values are enclosed in double-quotes (`"`) in this file. Double double-quotes are replaced with two double-quotes. The file contains the following fields:
movieId – movie id (49,174 movies)
title – movie title
releaseDate – the date when the movie was released
directedBy – directors separated by commas (`,`)
starring – cast separated by commas (`,`)
imdbId – id in IMDB
tmdbId – id in The movie DB
genres – genres separated by commas (`,`)

Tags Data File Structure (tags.csv)
-----------------------------------

All tags are contained in the file `tags.csv`. Each line of this file after the header row represents one tag applied to one movie by one user. Since some values contain commas (`,`), all the values are enclosed in double-quotes (`"`) in this file. Double double-quotes are replaced with two double-quotes. The contains the following fields (overall, there are 628,157 records):
userId - user id (8,908)
movieId – movie id (20,368)
tag - textual representation of the tag (45,198 tags)
timestamp – the tags were assigned to movies from January 22, 2007 till January 14, 2018

Tags are user-generated metadata about movies. Each tag is typically a single word or short phrase. The meaning, value and purpose of a particular tag is determined by each user.

Timestamps represent seconds since midnight Coordinated Universal Time (UTC) of January 1, 1970.

Tag Genome Data File Structure (tag_genome.csv)
-----------------------------------------

This file contains Tag Genome description for 11,005 movies from `answers.csv`, `recommendations.csv` and `training.csv`.

[genome-paper]: http://files.grouplens.org/papers/tag_genome.pdf

The tag genome is a data structure that contains tag relevance scores for movies.  The structure is a dense matrix: each movie in the genome has a value for *every* tag in the genome.

As described in [genome-paper], the tag genome encodes how strongly movies exhibit particular properties represented by tags (atmospheric, thought-provoking, realistic, etc.). The tag genome was computed using a machine learning algorithm on user-contributed content including tags, ratings, and textual reviews. The file contains 12,413,640 records and the following fields:
movieId – movie id (11,005 movies)
tag – textual representation of the tag (1,128 tags)
relevance – score which indicates how relevant the tag is to the movie on the scale from 0 to 1.

Please include the following citation if referencing tag genome data:

Jesse Vig, Shilad Sen, and John Riedl. 2012. The Tag Genome: Encoding Community Knowledge to Support Novel Interaction. ACM Trans. Interact. Intell. Syst. 2, 3: 13:1–13:44. https://doi.org/10.1145/2362394.2362395