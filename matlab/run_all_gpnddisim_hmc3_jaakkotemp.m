%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data/tempcodebranch/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/share/work/jtpelto/synergy-data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
%mybasedir_data='~/projects/pol2rnaseq/';
mybasedir_data='/share/synergy/';


mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/'];
% for lnDiffErfs.m
path4=[mybasedir_code 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir_code 'prior/matlab/'];
% for dist2.m
path6=[mybasedir_code 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir_code 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir_code 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)


load([mybasedir_data, 'data/series_for_matti_corrected_ver1.mat']);

changing_genes = importdata([mybasedir_data, 'analyses/changing_genes.txt']);


good_genes = { 'ENSG00000152413',
'ENSG00000175901',
'ENSG00000186628',
'ENSG00000217791',
'ENSG00000178401',
'ENSG00000230807',
'ENSG00000060303',
'ENSG00000235036',
'ENSG00000230207',
'ENSG00000249822',
'ENSG00000242013',
'ENSG00000247017',
'ENSG00000198042',
'ENSG00000152804',
'ENSG00000035403',
'ENSG00000198744',
'ENSG00000165863',
'ENSG00000105976',
'ENSG00000173542',
'ENSG00000135940',
'ENSG00000113369',
'ENSG00000233287',
'ENSG00000125962',
'ENSG00000120798',
'ENSG00000225920',
'ENSG00000165185',
'ENSG00000243404',
'ENSG00000133048',
'ENSG00000251536',
'ENSG00000128872',
'ENSG00000196678',
'ENSG00000213297',
'ENSG00000028310',
'ENSG00000234841',
'ENSG00000088325',
'ENSG00000134533',
'ENSG00000081665',
'ENSG00000019102',
'ENSG00000231610',
'ENSG00000154040',
'ENSG00000204308',
'ENSG00000152253',
'ENSG00000213073',
'ENSG00000073282',
'ENSG00000107562',
'ENSG00000251948',
'ENSG00000182718',
'ENSG00000186815',
'ENSG00000249660',
'ENSG00000140987',
'ENSG00000175602',
'ENSG00000125735',
'ENSG00000160813',
'ENSG00000179021',
'ENSG00000162852',
'ENSG00000249921',
'ENSG00000122965',
'ENSG00000164983',
'ENSG00000180769',
'ENSG00000138468',
'ENSG00000164463',
'ENSG00000105379',
'ENSG00000142178',
'ENSG00000139514',
'ENSG00000183808',
'ENSG00000226745',
'ENSG00000168300',
'ENSG00000181444',
'ENSG00000133639',
'ENSG00000173402',
'ENSG00000115306',
'ENSG00000116001',
'ENSG00000130649',
'ENSG00000168702',
'ENSG00000137992',
'ENSG00000175895',
'ENSG00000054965',
'ENSG00000163362',
'ENSG00000134897',
'ENSG00000113312',
'ENSG00000119321',
'ENSG00000183161',
'ENSG00000249969',
'ENSG00000187735',
'ENSG00000142207',
'ENSG00000068784',
'ENSG00000196843',
'ENSG00000196659',
'ENSG00000125652',
'ENSG00000052749',
'ENSG00000122779',
'ENSG00000085231',
'ENSG00000198691',
'ENSG00000213867',
'ENSG00000177728',
'ENSG00000138735',
'ENSG00000168101',
'ENSG00000241279',
'ENSG00000180953',
'ENSG00000163655' };

bad_inits = {'ENSG00000204652',
'ENSG00000248478',
'ENSG00000171711',
'ENSG00000232112',
'ENSG00000239195',
'ENSG00000248503',
'ENSG00000239208',
'ENSG00000239257',
'ENSG00000136149',
'ENSG00000186051',
'ENSG00000224871',
'ENSG00000232210',
'ENSG00000248547',
'ENSG00000186076',
'ENSG00000239283',
'ENSG00000085514',
'ENSG00000187140',
'ENSG00000205517',
'ENSG00000225125',
'ENSG00000232618',
'ENSG00000205527',
'ENSG00000248794',
'ENSG00000232703',
'ENSG00000225313',
'ENSG00000232713',
'ENSG00000160472',
'ENSG00000187266',
'ENSG00000232750',
'ENSG00000172638',
'ENSG00000205861',
'ENSG00000232751',
'ENSG00000239884',
'ENSG00000200879',
'ENSG00000230446',
'ENSG00000152439',
'ENSG00000230479',
'ENSG00000237173',
'ENSG00000201208',
'ENSG00000201302',
'ENSG00000222112',
'ENSG00000201376',
'ENSG00000222449',
'ENSG00000179168',
'ENSG00000250928',
'ENSG00000228559',
'ENSG00000250962',
'ENSG00000228561',
'ENSG00000235262',
'ENSG00000243955',
'ENSG00000146267',
'ENSG00000228598',
'ENSG00000251020',
'ENSG00000235290',
'ENSG00000243993',
'ENSG00000251063',
'ENSG00000202361',
'ENSG00000223669',
'ENSG00000230978',
'ENSG00000247421',
'ENSG00000223681',
'ENSG00000231028',
'ENSG00000237697',
'ENSG00000202487',
'ENSG00000231035',
'ENSG00000247508',
'ENSG00000203321',
'ENSG00000223741',
'ENSG00000237772',
'ENSG00000170409',
'ENSG00000203354',
'ENSG00000184502',
'ENSG00000231130',
'ENSG00000237782',
'ENSG00000216168',
'ENSG00000229211',
'ENSG00000229267',
'ENSG00000235732',
'ENSG00000244453',
'ENSG00000107796',
'ENSG00000216439',
'ENSG00000235796',
'ENSG00000244459',
'ENSG00000217159',
'ENSG00000229312',
'ENSG00000148082',
'ENSG00000244484',
'ENSG00000244496',
'ENSG00000251485',
'ENSG00000204377',
'ENSG00000224638',
'ENSG00000231951',
'ENSG00000238597',
'ENSG00000231966',
'ENSG00000238705',
'ENSG00000248412',
'ENSG00000231984',
'ENSG00000238820',
'ENSG00000232070',
'ENSG00000238881',
'ENSG00000232073',
'ENSG00000239017',
'ENSG00000248474',
'ENSG00000236376',
'ENSG00000252253',
'ENSG00000151131',
'ENSG00000199890',
'ENSG00000252311',
'ENSG00000200127',
'ENSG00000229947',
'ENSG00000245317',
'ENSG00000252494',
'ENSG00000200152',
'ENSG00000221164',
'ENSG00000236641',
'ENSG00000252531',
'ENSG00000200185',
'ENSG00000229970',
'ENSG00000252564',
'ENSG00000200259',
'ENSG00000218582',
'ENSG00000236071',
'ENSG00000199165',
'ENSG00000218631',
'ENSG00000229593',
'ENSG00000236076',
'ENSG00000219102',
'ENSG00000251790',
'ENSG00000167744',
'ENSG00000199347',
'ENSG00000229690',
'ENSG00000251791',
'ENSG00000199363',
'ENSG00000236151',
'ENSG00000251823',
'ENSG00000199366',
'ENSG00000219702',
'ENSG00000236209',
'ENSG00000217767',
'ENSG00000229380',
'ENSG00000148459',
'ENSG00000218073',
'ENSG00000229447',
'ENSG00000236001',
'ENSG00000244668',
'ENSG00000251603',
'ENSG00000218107',
'ENSG00000244723',
'ENSG00000251615',
'ENSG00000181123',
'ENSG00000218143',
'ENSG00000229487',
'ENSG00000251643',
'ENSG00000199053',
'ENSG00000218363',
'ENSG00000236060',
'ENSG00000251644',
'ENSG00000231587',
'ENSG00000185432',
'ENSG00000238118',
'ENSG00000248198',
'ENSG00000224219',
'ENSG00000248206',
'ENSG00000238137',
'ENSG00000231744',
'ENSG00000238168',
'ENSG00000248240',
'ENSG00000224251',
'ENSG00000238186',
'ENSG00000110169',
'ENSG00000230022',
'ENSG00000236677',
'ENSG00000200293',
'ENSG00000221319',
'ENSG00000230065',
'ENSG00000236703',
'ENSG00000252759',
'ENSG00000200394',
'ENSG00000221365',
'ENSG00000236732',
'ENSG00000245756',
'ENSG00000252865',
'ENSG00000200428',
'ENSG00000221410',
'ENSG00000252926',
'ENSG00000221461',
'ENSG00000253047',
'ENSG00000200502',
'ENSG00000230098',
'ENSG00000253096',
'ENSG00000199382',
'ENSG00000229747',
'ENSG00000199426',
'ENSG00000236263',
'ENSG00000251992',
'ENSG00000181735',
'ENSG00000199645',
'ENSG00000236307',
'ENSG00000252174',
'ENSG00000199649',
'ENSG00000220014',
'ENSG00000252192',
'ENSG00000199697',
'ENSG00000220305',
'ENSG00000229853',
'ENSG00000236325',
'ENSG00000182057',
'ENSG00000199719',
'ENSG00000229862',
'ENSG00000252207',
'ENSG00000166667',
'ENSG00000235573',
'ENSG00000251234',
'ENSG00000229052',
'ENSG00000215856',
'ENSG00000251254',
'ENSG00000129455',
'ENSG00000215871',
'ENSG00000235629',
'ENSG00000244307',
'ENSG00000216031',
'ENSG00000235672',
'ENSG00000244349',
'ENSG00000229190',
'ENSG00000235681',
'ENSG00000244416',
'ENSG00000203417',
'ENSG00000223813',
'ENSG00000231177',
'ENSG00000203437',
'ENSG00000237857',
'ENSG00000247638',
'ENSG00000184825',
'ENSG00000203486',
'ENSG00000223849',
'ENSG00000237892',
'ENSG00000247694',
'ENSG00000203542',
'ENSG00000223861',
'ENSG00000231327',
'ENSG00000237911',
'ENSG00000203564',
'ENSG00000231340',
'ENSG00000200579',
'ENSG00000221500',
'ENSG00000230122',
'ENSG00000200591',
'ENSG00000221526',
'ENSG00000230126',
'ENSG00000111291',
'ENSG00000182873',
'ENSG00000221667',
'ENSG00000230171',
'ENSG00000237015',
'ENSG00000200795',
'ENSG00000221879',
'ENSG00000230200',
'ENSG00000237039',
'ENSG00000049247',
'ENSG00000111348',
'ENSG00000200851',
'ENSG00000099365',
'ENSG00000121897',
'ENSG00000207696',
'ENSG00000226104',
'ENSG00000240808',
'ENSG00000188765',
'ENSG00000207733',
'ENSG00000240825',
'ENSG00000207741',
'ENSG00000233476',
'ENSG00000207783',
'ENSG00000226176',
'ENSG00000233538',
'ENSG00000207790',
'ENSG00000240996',
'ENSG00000249412',
'ENSG00000196476',
'ENSG00000233996',
'ENSG00000226920',
'ENSG00000226945',
'ENSG00000234070',
'ENSG00000241821',
'ENSG00000249821',
'ENSG00000140749',
'ENSG00000213082',
'ENSG00000226957',
'ENSG00000234261',
'ENSG00000241907',
'ENSG00000249855',
'ENSG00000176933',
'ENSG00000197049',
'ENSG00000234473',
'ENSG00000250068',
'ENSG00000142227',
'ENSG00000197061',
'ENSG00000234498',
'ENSG00000242567',
'ENSG00000250144',
'ENSG00000124939',
'ENSG00000227368',
'ENSG00000250162',
'ENSG00000213403',
'ENSG00000227407',
'ENSG00000234535',
'ENSG00000242639',
'ENSG00000071677',
'ENSG00000185631',
'ENSG00000224261',
'ENSG00000231795',
'ENSG00000238240',
'ENSG00000248253',
'ENSG00000224314',
'ENSG00000231835',
'ENSG00000238279',
'ENSG00000135547',
'ENSG00000185672',
'ENSG00000224450',
'ENSG00000231856',
'ENSG00000231888',
'ENSG00000238363',
'ENSG00000224468',
'ENSG00000238387',
'ENSG00000248388',
'ENSG00000201405',
'ENSG00000237259',
'ENSG00000201678',
'ENSG00000222493',
'ENSG00000230642',
'ENSG00000223361',
'ENSG00000230679',
'ENSG00000112149',
'ENSG00000201772',
'ENSG00000183598',
'ENSG00000223401',
'ENSG00000237434',
'ENSG00000201823',
'ENSG00000223414',
'ENSG00000201838',
'ENSG00000230834',
'ENSG00000202058',
'ENSG00000223451',
'ENSG00000247313',
'ENSG00000202077',
'ENSG00000237520',
'ENSG00000202186',
'ENSG00000230964',
'ENSG00000237575',
'ENSG00000202314',
'ENSG00000230970',
'ENSG00000237641',
'ENSG00000202360',
'ENSG00000230974',
'ENSG00000237659',
'ENSG00000228873',
'ENSG00000235419',
'ENSG00000198671',
'ENSG00000215791',
'ENSG00000235448',
'ENSG00000198637',
'ENSG00000235423',
'ENSG00000251206',
'ENSG00000232470',
'ENSG00000239607',
'ENSG00000224967',
'ENSG00000232482',
'ENSG00000239617',
'ENSG00000239688',
'ENSG00000225094',
'ENSG00000232578',
'ENSG00000239718',
'ENSG00000225113',
'ENSG00000232608',
'ENSG00000239726',
'ENSG00000203606',
'ENSG00000231365',
'ENSG00000237969',
'ENSG00000185065',
'ENSG00000223998',
'ENSG00000231373',
'ENSG00000224014',
'ENSG00000231381',
'ENSG00000237979',
'ENSG00000203814',
'ENSG00000238015',
'ENSG00000248157',
'ENSG00000204044',
'ENSG00000224080',
'ENSG00000231581',
'ENSG00000171943',
'ENSG00000239302',
'ENSG00000204894',
'ENSG00000232334',
'ENSG00000239367',
'ENSG00000248577',
'ENSG00000232354',
'ENSG00000239453',
'ENSG00000224911',
'ENSG00000239483',
'ENSG00000248590',
'ENSG00000136327',
'ENSG00000204982',
'ENSG00000224934',
'ENSG00000239516',
'ENSG00000172186',
'ENSG00000248639',
'ENSG00000207501',
'ENSG00000233148',
'ENSG00000207523',
'ENSG00000225969',
'ENSG00000233197',
'ENSG00000240667',
'ENSG00000173838',
'ENSG00000207547',
'ENSG00000225981',
'ENSG00000233202',
'ENSG00000249292',
'ENSG00000207611',
'ENSG00000225991',
'ENSG00000240690',
'ENSG00000249313',
'ENSG00000207617',
'ENSG00000233337',
'ENSG00000207623',
'ENSG00000226080',
'ENSG00000233391',
'ENSG00000161570',
'ENSG00000187959',
'ENSG00000206874',
'ENSG00000225770',
'ENSG00000232986',
'ENSG00000137875',
'ENSG00000206941',
'ENSG00000240322',
'ENSG00000225790',
'ENSG00000207112',
'ENSG00000225803',
'ENSG00000240366',
'ENSG00000207129',
'ENSG00000225811',
'ENSG00000233021',
'ENSG00000207145',
'ENSG00000233025',
'ENSG00000225476',
'ENSG00000240125',
'ENSG00000249055',
'ENSG00000206634',
'ENSG00000225518',
'ENSG00000232921',
'ENSG00000137726',
'ENSG00000206728',
'ENSG00000249076',
'ENSG00000206834',
'ENSG00000232937',
'ENSG00000240206',
'ENSG00000249087',
'ENSG00000161328',
'ENSG00000206838',
'ENSG00000232969',
'ENSG00000240232',
'ENSG00000206861',
'ENSG00000240270',
'ENSG00000249121',
'ENSG00000228638',
'ENSG00000251073',
'ENSG00000166268',
'ENSG00000228643',
'ENSG00000235313',
'ENSG00000244063',
'ENSG00000166321',
'ENSG00000179467',
'ENSG00000215462',
'ENSG00000228661',
'ENSG00000244078',
'ENSG00000251161',
'ENSG00000187475',
'ENSG00000225379',
'ENSG00000239912',
'ENSG00000160886',
'ENSG00000225450',
'ENSG00000232824',
'ENSG00000248993',
'ENSG00000206573',
'ENSG00000187624',
'ENSG00000206622',
'ENSG00000240050',
'ENSG00000249044',
'ENSG00000228741',
'ENSG00000235343',
'ENSG00000244094',
'ENSG00000251172',
'ENSG00000228857',
'ENSG00000235369',
'ENSG00000244121',
'ENSG00000251174',
'ENSG00000207166',
'ENSG00000225842',
'ENSG00000207291',
'ENSG00000240441',
'ENSG00000207311',
'ENSG00000225913',
'ENSG00000173588',
'ENSG00000207314',
'ENSG00000233069',
'ENSG00000249251',
'ENSG00000207395',
'ENSG00000225941',
'ENSG00000233108',
'ENSG00000240527',
'ENSG00000173727',
'ENSG00000225963',
'ENSG00000212588',
'ENSG00000233893',
'ENSG00000241434',
'ENSG00000249681',
'ENSG00000140092',
'ENSG00000233899',
'ENSG00000249728',
'ENSG00000175906',
'ENSG00000233903',
'ENSG00000249741',
'ENSG00000140398',
'ENSG00000226823',
'ENSG00000241634',
'ENSG00000207948',
'ENSG00000241157',
'ENSG00000100079',
'ENSG00000174912',
'ENSG00000207981',
'ENSG00000226276',
'ENSG00000241170',
'ENSG00000208028',
'ENSG00000226312',
'ENSG00000241202',
'ENSG00000249526',
'ENSG00000189385',
'ENSG00000226424',
'ENSG00000241247',
'ENSG00000249531',
'ENSG00000163202',
'ENSG00000209707',
'ENSG00000226457',
'ENSG00000226507',
'ENSG00000249548',
'ENSG00000100271',
'ENSG00000211513',
'ENSG00000233799',
'ENSG00000139800',
'ENSG00000211553',
'ENSG00000233818',
'ENSG00000241361',
'ENSG00000196110',
'ENSG00000211591',
'ENSG00000241362',
'ENSG00000249572',
'ENSG00000241413',
'ENSG00000249615',
'ENSG00000212296',
'ENSG00000233871',
'ENSG00000226521',
'ENSG00000241358',
'ENSG00000234353',
'ENSG00000249996',
'ENSG00000242154',
'ENSG00000250030',
'ENSG00000124693',
'ENSG00000227268',
'ENSG00000234450',
'ENSG00000227297',
'ENSG00000234465',
'ENSG00000235063',
'ENSG00000243597',
'ENSG00000197882',
'ENSG00000243629',
'ENSG00000235066',
'ENSG00000228350',
'ENSG00000250731',
'ENSG00000214894',
'ENSG00000177236',
'ENSG00000197084',
'ENSG00000213412',
'ENSG00000227440',
'ENSG00000242675',
'ENSG00000227452',
'ENSG00000234552',
'ENSG00000242714',
'ENSG00000250229',
'ENSG00000227505',
'ENSG00000250255',
'ENSG00000213558',
'ENSG00000227541',
'ENSG00000242756',
'ENSG00000250258',
'ENSG00000243179',
'ENSG00000250535',
'ENSG00000126264',
'ENSG00000126266',
'ENSG00000178550',
'ENSG00000197635',
'ENSG00000214212',
'ENSG00000234975',
'ENSG00000213657',
'ENSG00000234630',
'ENSG00000242833',
'ENSG00000250301',
'ENSG00000242959',
'ENSG00000250321',
'ENSG00000227666',
'ENSG00000242992',
'ENSG00000250324',
'ENSG00000227694',
'ENSG00000234724',
'ENSG00000228415',
'ENSG00000243738',
'ENSG00000228439',
'ENSG00000250796',
'ENSG00000228487',
'ENSG00000235167',
'ENSG00000105538',
'ENSG00000165935',
'ENSG00000243015',
'ENSG00000250410',
'ENSG00000227741',
'ENSG00000234775',
'ENSG00000243054',
'ENSG00000250419',
'ENSG00000125787',
'ENSG00000227896',
'ENSG00000234776',
'ENSG00000243055',
'ENSG00000250474',
'ENSG00000243071',
'ENSG00000197522',
'ENSG00000176236',
'ENSG00000241954',
'ENSG00000241984',
'ENSG00000249911',
'ENSG00000234311',
'ENSG00000241988',
'ENSG00000124455',
'ENSG00000141294',
'ENSG00000234315',
'ENSG00000101440',
'ENSG00000227192',
'ENSG00000176472',
'ENSG00000234345',
'ENSG00000214239',
'ENSG00000228124',
'ENSG00000234981',
'ENSG00000143942',
'ENSG00000197723',
'ENSG00000243352',
'ENSG00000243431',
'ENSG00000228216',
'ENSG00000243483',
'ENSG00000235028',
'ENSG00000250679' };

good_genes = bad_inits;



mygenes = changing_genes;
indices = zeros(size(mygenes));  % index of the desired genes in the bininfo array
for k=1:length(mygenes),
  t = mygenes{k};
  I = find(str2double(t(10:end)) == bininfo(:, 5));
  if ~isempty(I),
    indices(k) = I;
  else
    indices(k) = NaN;
  end
end





%---------------------------------------------------
% Load maximum likelihood fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)

%load allresults_shifted_polorrna4.mat
% allresults_ensemblids_pol2 allresults_geneindices_pol2 allresults_jointmodels_pol2 allresults_jointtransforminfos_pol2 allresults_loglikelihoods_pol2
% allresults_ensemblids_rna allresults_geneindices_rna allresults_jointmodels_rna allresults_jointtransforminfos_rna allresults_loglikelihoods_rna
%allresults_jointtransforminfos_pol2 = allresults_tinfos_pol2;
%allresults_jointtransforminfos_rna = allresults_tinfos_rna;


%load allresults_shifted_longerdelay5.mat
%allresults_ensemblids_joint = allresults_ensemblids; 
%allresults_geneindices_joint = allresults_geneindices;
%allresults_jointmodels_joint = allresults_jointmodels;
%allresults_jointtransforminfos_joint = allresults_jointtransforminfos;
%allresults_loglikelihoods_joint = allresults_loglikelihoods;




%---------------------------------------------------
% Compute estimated probability of joint fit being useful
%---------------------------------------------------

% probcomparison = nan*ones(max([length(allresults_jointtransforminfos_pol2) ...
%      length(allresults_jointtransforminfos_joint)]),1);
% for i=1:length(probcomparison),
%   if (i<=size(allresults_loglikelihoods_joint,1)) ...
%     && (i<=size(allresults_loglikelihoods_pol2,1)),
%     if (~isempty(allresults_jointmodels_joint{i})) && (~isempty(allresults_jointmodels_pol2{i})),
%       probcomparison(i)=allresults_loglikelihoods_joint(i,3) ...
%           -allresults_loglikelihoods_pol2(i,3) -allresults_loglikelihoods_rna(i,3);
%     end;
%   end;  
% end;


% Rank genes by probability of joint fit being useful
%[y,Ijointrank]=sort(-probcomparison);


%maxHMCgenes = 10;
nHMCiters = 10000;
%HMCsamples = cell(length(allresults_jointmodels_joint),1);

%istart=1;
%iend=maxHMCgenes;
myI=1:length(indices);

for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  gene_index = indices(i);
  gene_name = mygenes{i};

  if isnan(gene_index),
    continue;
  end

  fprintf('Running gene %d/%d: %s\n', find(i==myI), length(myI), gene_name);

  randn('seed',bininfo(gene_index,5));
  rand('seed',bininfo(gene_index,5)+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  lengthscale=2;

  timevector = {timevector,timevector};
  dataVals = {tempvals1, tempvals2};
  initializationtype=1;

  try,
    [m,temptransforminfo]=createNdSimDisim_celltimes(timevector, ...
						     dataVals,lengthscale,initializationtype,[],[],1);
  catch,
    warning('Unable to create a model for gene %s\n', gene_name);
    continue;
  end

  % Apply HMC sampling
  HMCsamples = gpnddisimSampleHMC(m, 1, nHMCiters);
  HMCsamples = HMCsamples(10:10:end, :);
  save(sprintf('hmc_results/%s_samples_%s.mat', gene_name, id), ...
       'gene_name', 'gene_index', 'm', 'HMCsamples');
end;

