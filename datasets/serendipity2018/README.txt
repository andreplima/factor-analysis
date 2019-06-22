#------------------------------------------------------------------------------------------------------------------------------
# This document comprises segments from the original README file of the Serendipity'2018 dataset, which has
# generously been made available by the authors:
# Kotkov, D., Konstan, J. A., Zhao, Q., & Veijalainen, J. (2018, April). Investigating serendipity in
#   recommender systems based on real user feedback. In Proceedings of the 33rd Annual ACM Symposium on
#   Applied Computing (pp. 1341-1350). ACM.
# The complete dataset is available at: https://grouplens.org/datasets/serendipity-2018/
#------------------------------------------------------------------------------------------------------------------------------

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

