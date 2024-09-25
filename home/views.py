from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from .forms import FastaUploadForm
from .services import PROSITEProcessingService, ListProcessingService, FastaService # Certifique-se de que o serviço está importado
from .aminoacid_colors import AminoacidColorMap

# Criar uma instância do seu serviço
fasta_service = FastaService()
list_processing_service = ListProcessingService()
prosite_processing_service = PROSITEProcessingService(fasta_service=fasta_service,list_processing_service=list_processing_service)

def home(request):
    title = "BIOINFORMÁTICA ESTRUTURAL"
    form = FastaUploadForm()
    context = {
        'title': title,
        'form': form,
    }
    return render(request, 'home/home.html', context)

@csrf_exempt
def upload_fasta(request):
    if request.method == 'POST':
        uploaded_file = request.FILES['fasta_file']  # Nome do campo do formulário
        contents = uploaded_file.read().decode('utf-8')
        
        # Obter valores do formulário
        score_model_conservation = request.POST.get('score_model_conservation')  # Obtém o valor do modelo de conservação
        xthreshold = request.POST.get('xthreshold')  # Obtém o valor do X-Threshold
        xthreshold = int(xthreshold) if xthreshold else None  # Converte para inteiro, se aplicável

        # Processar o conteúdo FASTA
        fasta_entries = prosite_processing_service.parseFasta(contents)
        prosite_signatures = prosite_processing_service.processFasta(
            contents,
            score_model_conservation,
            xthreshold
        )

        return JsonResponse({
            'fasta_entries': fasta_entries,
            'prosite_signatures': prosite_signatures,
        })
    return JsonResponse({'error': 'Método não permitido.'}, status=405)

def get_color(aminoacid):
    return AminoacidColorMap.COLOR_MAP_HEX[aminoacid]

def fasta_view(request):
    fasta_entries = []
    prosite_assinatures = []
    if request.method == 'POST' and request.FILES['fasta_file']:
        file = request.FILES['fasta_file']
        contents = file.read().decode('utf-8')
        
        prosite_service = PROSITEProcessingService()
        fasta_entries = prosite_service.parse_fasta(contents)
        prosite_assinatures = prosite_service.process_fasta(contents, "some_model", 20)  # Ajuste conforme necessário

    color_map_hex = AminoacidColorMap.COLOR_MAP_HEX
    context = {
        'fasta_entries': fasta_entries,
        'prosite_assinatures': prosite_assinatures,
        'color_map': color_map_hex,
    }
    return render(request, 'fasta_display.html', context)