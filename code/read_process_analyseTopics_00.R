######################### Load and analyze bibtex file containing the publications by AJMS ###################
#################### directory structure #############################################
cDir = '/dyStorage/Analyze_topics_ajms/codes'
dDir = '/dyStorage/Analyze_topics_ajms/rawData'
docsDir = '/dyStorage/Analyze_topics_ajms/docs'
outFDir = '/dyStorage/Analyze_topics_ajms/outFig_tabs'

############################### load required libraries #########################
lib_req = c('tidytext','tidyr','dplyr','ggplot2','bib2df','stringr','ldatuning','topicmodels')
lapply(lib_req,library,character.only=T)

############################# Load the raw data ###################################
setwd(dDir)

bib_ajms = bib2df('Afr_jMar_pubs.bib')

## create better journal code
bib_ajms$new_journal = ifelse(str_detect(string = str_to_lower(bib_ajms$JOURNAL),pattern = '^afr'),'AJMS',
                              ifelse(str_detect(string = str_to_lower(bib_ajms$JOURNAL),pattern = '^sou'),'SAJMS',NA))

## select relevant variables and unnest
bib_ajms %>% 
  select(YEAR, AUTHOR,TITLE,PAGES,ABSTRACT)%>% 
  unnest()->bib_ajmsUnNest

## remove those without journal information, and missing abstract
bib_ajms%>%filter(!is.na(bib_ajms$new_journal))%>%select(YEAR, ABSTRACT)%>%rename(Year=YEAR,Abstract=ABSTRACT)%>%
  arrange(desc(-Year))%>%filter(!is.na(Abstract))->bib_ajmsSel
tab_yr =  table(bib_ajmsSel$Year)
nP_yrs = sapply(1:length(unique(bib_ajmsSel$Year)), function(x,nP) {paste('p',1:nP[x],sep='')},nP=tab_yr)%>%unlist()
bib_ajmsSel<-bib_ajmsSel%>%mutate(Paper = nP_yrs)%>%select(Year,Paper,Abstract)

####### Check the total number of publications by year
bib_nPprs = bib_ajmsSel%>%group_by(Year)%>%count(Year)%>%select(Year,n)

##################### bigrams ######################
#### prepare additional stop words to remove
add_stop = data.frame(word=c('en','cent','mm','m','al','textperiodcentered','van','de',
                             'cm','nisc',"pty",'sd','ac',"lg","backslash",'textless','textgreater','textcopyright'),
                      stringsAsFactors = F)

## create bigrams (pair of word that occur in sequence)
bib_ajmsTidy_bigram = bib_ajmsSel%>%unnest_tokens(bigram,Abstract,token = 'ngrams',n=2)

## prepare the bigrams for cleaning (removing stopwords)
bib_ajmsTidy_bigram2 <- bib_ajmsTidy_bigram %>%
  separate(bigram, c("word1", "word2"), sep = " ")
bigrams_filtered <- bib_ajmsTidy_bigram2 %>%
  filter(!word1 %in% stop_words$word) %>%
  filter(!word2 %in% stop_words$word) %>%filter(!word1 %in% add_stop$word)%>%
  filter(!word2 %in% add_stop$word) %>%filter(!str_detect(word1,'[[:digit:]]|[[:punct:]]'))%>%
  filter(!str_detect(word2,'[[:digit:]]|[[:punct:]]'))


# new bigram counts:
bigram_counts <- bigrams_filtered %>%
  count(word1, word2, sort = TRUE)

## combining the final sets of bigrams
bigrams_united <- bigrams_filtered %>%
  unite(bigram, word1, word2, sep = " ")
bigrams_united

bigrams_unitedCopunt <- bigrams_united%>%count(bigram,sort = TRUE)

############## tidytext processing of the abstracts 
bib_ajmsSel$document <- paste(bib_ajmsSel$Year,bib_ajmsSel$Paper,sep = '_')
bib_ajmsTidy = bib_ajmsSel%>%unnest_tokens(word,Abstract)

########### clean text (removing stop words)
bib_clean1 <- bib_ajmsTidy%>%anti_join(stop_words)
bib_clean1a <- bib_clean1%>%anti_join(add_stop);bib_clean1a$word[bib_clean1a$word=='demersus']='demersal'
bib_clean1b <-bib_clean1a[!str_detect(bib_clean1a$word,'[[:digit:]]|[[:punct:]]'),]
bib_clean2 <- bib_clean1b[,-c(1,2)]%>%group_by(document)%>%filter(word!='')%>%count(word,sort = TRUE)

bib_ajmsSel_2 = bib_ajmsSel[bib_ajmsSel$document%in%unique(bib_clean2$document),]
# convert the cleaned text document to document term matrix 
bib_dtm <-bib_clean2%>%cast_dtm(document,word,n)

###  word count for the overall word cloud 
bib_cleanAll <- bib_clean1b[,-c(1,2)]%>%filter(word!='')%>%count(word,sort = TRUE)

##### Find optmal number of topics 

n_tops = 5:50
Mtrics = "CaoJuan2009" 
system.time(opt_topics <- FindTopicsNumber(dtm = bib_dtm,topics = n_tops,metrics = Mtrics,mc.cores = 6,verbose = TRUE))

### Plot profile of the performance measure 
tiff('optimal_no_topics_caoJuan.tiff',units = 'cm',height = 25,width = 15,res = 300,compression = 'lzw')
FindTopicsNumber_plot(opt_topics[,1:2])+
  theme(strip.text.x = element_text(size=6),axis.text = element_text(size=5))
dev.off()

### LDA model using the optimal number of K identified above
system.time(top_24 <- LDA(x = bib_dtm,k = 24,control = list(seed = 1235)))

################# Save model selection result and the final model fitted ####################

save(list=c('opt_topics','top_24'),file='model_selectionFit.RData')

################## tidy and prepare for visualization #######################
setwd(outFDir)
load('model_selectionFit.RData')

### word probability for each topic and term combination (probability a word/term being generated from a topic) (beta)
tidy_top_24b <-tidy(top_24,matrix='beta')

## As LDA models each document as being made of mixture of topics and each topic being made of a mixture of terms/word.
## We can also compute the probability for each combination of document and topic. The probability that terms/words in
## each document are from a particular topic (gamma). Can be seen each document is made of mixture different proportion 
## of the different topics

tidy_top_24g <- tidy(top_24,matrix='gamma')

## processing extracted results (gamma)

abs_top24_gamma <-tidy_top_24g%>%
  separate(document,c("Year","paper"),sep = '_',convert = TRUE)

## processing extracted results (beta)
abs_top24_beta <- tidy_top_24b %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic,-beta)
abs_top24_beta


## extract the sets of words that characterize each of the topices (in this case just five words)
terms24 = terms(top_24,5)
topics24 = topics(top_24)

## posterior probability of topics in each document and terms for each topics
posterior_24 <- posterior(top_24) 

wt_24 <- as.data.frame(t(posterior_24$terms)) # Posterior probabilitY of terms for each topics
wt_24L <- wt_24 %>%
  mutate(word = rownames(wt_24)) %>%
  tidyr::gather(topic, weight, -word) 

######## Frequency of occurrence of topics  
terms_24 <- apply(terms24, MARGIN = 2, paste, collapse = ", ")

## merge each document index to the topics assigned based on the LDA model fitted above
topicsAbs_24 <- data.frame(document=bib_ajmsSel_2$document, topic=topics24)

### a bit of processing split the document (unique index) into a year and paper 
topicAbs_24_sep <- topicsAbs_24%>%separate(document,c('Year','paper'),sep = '_',convert = TRUE)

### top five terms that characterize each of the 24 topics
terms24_df = terms_24%>%as.data.frame(stringsAsFactors=FALSE)%>%rename(keys = '.')%>%
  mutate(topic = 1:24,topic_key=paste(topic,": ",keys,sep=''))
## Topics into which each of the documents belong to and the corresponding terms that characterize each topics
topicAbs_24_sep_j = topicAbs_24_sep%>%left_join(terms24_df[,c(2:3)])%>%arrange(desc(topic))
topicAbs_24_sep_j$topic_key <-reorder(factor(topicAbs_24_sep_j$topic_key),X = topicAbs_24_sep_j$topic)


############# some visual summary of the time series of the 24  topics 
ggplot(topicAbs_24_sep,aes(Year,..count..))+geom_density(aes(fill=terms_24[topic]),alpha=0.5,position='identity')+theme_bw()

###### slightly different variants of same figures as above
ggplot(topicAbs_24_sep,aes(Year,..count..))+geom_bar()+facet_wrap(~terms_24[topic],nrow = 5)+theme_bw()+
  theme(strip.text.x = element_text(size=6))-> topic_byYear1
ggplot(topicAbs_24_sep_j,aes(Year,..count..))+geom_bar()+facet_wrap(~topic_key,nrow = 5)+theme_bw()+
  theme(strip.text.x = element_text(size=5.5))-> topic_byYear1a
ggplot(topicAbs_24_sep,aes(Year,..count..))+geom_bar()+facet_wrap(~topic,nrow = 5)+theme_bw() -> topic_byYear2


ggsave(filename = 'time_series_topics_ver1.tiff',plot = topic_byYear1,device = 'tiff',width = 25,height = 15,units = 'cm',
       dpi = 300,compression='lzw')
ggsave(filename = 'time_series_topics_ver1a.tiff',plot = topic_byYear1a,device = 'tiff',width = 25,height = 15,units = 'cm',
       dpi = 300,compression='lzw')
ggsave(filename = 'time_series_topics_ver2.tiff',plot = topic_byYear2,device = 'tiff',width = 25,height = 15,units = 'cm',
       dpi = 300,compression='lzw')

#### plots for beta ()
abs_top24_beta %>%
  mutate(term = reorder(term, beta)) %>%
  ggplot(aes(term, beta, fill = factor(topic))) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  facet_wrap(~ topic, scales = "free") +
  coord_flip()

############# word clouds by topic ####################
setwd(outFDir)

library(animation)
library(wordcloud)
library(RColorBrewer)

n=100

palette='Greens'
pal <- rep(brewer.pal(9, palette), each = ceiling(n / 9))[n:1]


par(mar=c(1,1,3,1),mgp=c(0.2,0.2,1))
saveGIF({for(i in 1:24){
  #tiff('word_clouds_for_topic_',i)
  wt3_L <- wt_24L %>%
    filter(topic == i) %>%
    arrange(desc(weight))
  with(wt3_L[1:n, ], 
       wordcloud(word, freq = weight, random.order = FALSE, 
                 ordered.colors = TRUE, colors = pal))
  mtext(paste('Topic_',i,sep=''), cex = 0.6,col = 'blue')  
}},'Optimally_identified_topics_AJMS.gif',interval=2)

dev.off()

n2=50
palette='Greens'
pal2 <- rep(brewer.pal(9, palette), each = ceiling(n2 / 9))[n2:1]

tiff('Word_clouds_24_topics_ajms.tif',unit='cm',height=18,width = 20,res = 350,compression = 'lzw')
par(mfrow=c(4,6),mar=c(0,1,3,1),mgp=c(0.2,0.2,1))
for(i in 1:24){
  #tiff('word_clouds_for_topic_',i)
  wt3_L <- wt_24L %>%
    filter(topic == i) %>%
    arrange(desc(weight))
  with(wt3_L[1:n2, ], 
       wordcloud(word, freq = weight, random.order = FALSE, scale = c(2,0.2),
                 ordered.colors = TRUE, colors = pal2))
  mtext(paste('Topic_',i,sep=''), cex = 0.6,col = 'blue')  
}
dev.off()

#### same as above but split into two
n2a=100
palette='Greens'
pal2a <- rep(brewer.pal(9, palette), each = ceiling(n2a / 9))[n2a:1]

par(mfrow=c(3,4),mar=c(0,1,3,1),mgp=c(0.2,0.2,1))
for(i in 1:12){
  #tiff('word_clouds_for_topic_',i)
  wt3_L <- wt_24L %>%
    filter(topic == i) %>%
    arrange(desc(weight))
  with(wt3_L[1:n2a, ], 
       wordcloud(word, freq = weight, random.order = FALSE, scale = c(2,0.2),
                 ordered.colors = TRUE, colors = pal2a))
  mtext(paste('Topic_',i,sep=''), cex = 0.6,col = 'blue')  
}


par(mfrow=c(3,4),mar=c(0,1,3,1),mgp=c(0.2,0.2,1))
for(i in 13:24){
  #tiff('word_clouds_for_topic_',i)
  wt3_L <- wt_24L %>%
    filter(topic == i) %>%
    arrange(desc(weight))
  with(wt3_L[1:n2a, ], 
       wordcloud(word, freq = weight, random.order = FALSE, scale = c(2,0.2),
                 ordered.colors = TRUE, colors = pal2a))
  mtext(paste('Topic_',i,sep=''), cex = 0.6,col = 'blue')  
}

################ word cloud based on all abstracts ####################
n3=400
pal3 <- rep(brewer.pal(9, 'Greens'), each = ceiling(n3 / 9))[n3:1]

tiff('Word_clouds_all_abstract_ajms.tif',unit='cm',height=15,width = 15,res = 350,compression = 'lzw')
par(mar=c(1,1,3,1),mgp=c(0.2,0.2,1))
with(bib_cleanAll[1:n3, ], 
     wordcloud(word, freq = n, random.order = FALSE, 
               ordered.colors = TRUE,colors = pal3))

dev.off()

################ word colouds based on bigrams #########################

n3=400
pal3 <- rep(brewer.pal(9, 'Greens'), each = ceiling(n3 / 9))[n3:1]

setwd(outFDir)
tiff('bigram_clouds_all_abstract_ajms.tif',unit='cm',height=15,width = 15,res = 350,compression = 'lzw')
par(mar=c(1,1,3,1),mgp=c(0.2,0.2,1))
with(bigrams_unitedCopunt[1:n3, ], 
     wordcloud(bigram, freq = n, random.order = FALSE, 
               ordered.colors = TRUE,colors = pal3))

dev.off()
